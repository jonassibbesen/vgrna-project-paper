
# plot_variant_concordance_effect.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/variant_r2/")

########

parse_exp_data <- function(filename) {
  
  print(filename)
  name_split <- strsplit(basename(filename), "_")[[1]]
  
  stopifnot(grepl("ENCSR", name_split[6]))
  
  data <- read_table(filename, col_names = T) %>%
    mutate(Chrom = as.character(Chrom))

  data <- data %>%    
    group_by(Chrom, Position) %>%
    mutate(var_exp = sum(TPM)) %>%
    filter(var_exp > 0) %>%
    add_column(Sample = name_split[6])
  
  return(data)
}


calc_count_inner <- function(min_prob, min_exp, exp_data) {
  
  print(min_prob)
  print(min_exp)
  
  exp_count <- exp_data %>% 
    filter(Probability >= min_prob & TPM >= min_exp) %>%
    group_by(Chrom, Position, AlleleNum, AlleleType, HomopolymerLength) %>%
    summarise(num_exp = n()) %>%
    mutate(AlleleType = as.factor(AlleleType)) %>%
    add_column(min_prob = min_prob) %>%
    add_column(min_exp = min_exp)
  
  return(exp_count)
}

calc_cons_inner <- function(min_prob, min_exp, exp_data) {
  
  print(min_prob)
  print(min_exp)
  
  exp_cons <- exp_data %>% 
    filter(var_exp >= min_exp) %>%
    mutate(is_exp = (Probability >= min_prob & TPM > 0)) %>%
    group_by(Chrom, Position, AlleleNum, AlleleType, HomopolymerLength) %>%
    count(is_exp) %>%
    rename(count = n) %>%
    mutate(num_count = n(), max_count = max(count), sum_count = sum(count)) %>%
    mutate(AlleleType = as.factor(AlleleType)) %>%
    add_column(min_prob = min_prob) %>%
    add_column(min_exp = min_exp)
  
  return(exp_cons)
}

calc_count <- function(min_exp, exp_data) {
  
  exp_count <- map_dfr(c(0.8), calc_count_inner, min_exp, exp_data)
  
  return(exp_count)
}

calc_cons <- function(min_exp, exp_data) {
  
  exp_cons <- map_dfr(c(0.8), calc_cons_inner, min_exp, exp_data)
  
  return(exp_cons)
}

parse_and_analyse_exp_data <- function(chrom) {
  
  gc()
  
  print(chrom)  
  
  exp_data <- map_dfr(list.files(path = paste("allele_expression/", chrom, sep = ""), pattern = "allele_exp_rpvg_mpmap.*.txt", full.names = T, recursive = T), parse_exp_data)
  
  print(unique(exp_data$Sample))
  
  if (chrom != "Y") {

    exp_count <- map_dfr(c(0.01, 1, 2, 5, 10, 15, 20), calc_count, exp_data)
    exp_cons <- map_dfr(c(0.01, 1, 2, 5, 10, 15, 20), calc_cons, exp_data)
    
  } else {
    
    exp_count <- map_dfr(c(0.01), calc_count, exp_data)
    exp_cons <- map_dfr(c(0.01), calc_cons, exp_data)
  }
  
  exp_count <- exp_count %>%
    group_by(AlleleType, HomopolymerLength, num_exp, min_prob, min_exp) %>%
    summarise(n = n())
  
  exp_cons <- exp_cons %>%
    group_by(AlleleType, HomopolymerLength, is_exp, count, num_count, sum_count, max_count, min_prob, min_exp) %>%
    summarise(n = n())
              
  save(exp_count, file = paste("rdata/exp_count_", chrom, ".RData", sep = "", collapse = ""))
  save(exp_cons, file = paste("rdata/exp_cons_", chrom, ".RData", sep = "", collapse = ""))
}

map_dfr(c(seq(1, 22), "X", "Y"), parse_and_analyse_exp_data)

exp_count_all <- list()

for (f in list.files(pattern = "exp_count_.*RData", full.names = T, recursive = T)) { 
  
  print(f)
  load(f)
  
  exp_count_all[[f]] <- exp_count
}

exp_cons_all <- list()

for (f in list.files(pattern = "exp_cons_.*RData", full.names = T, recursive = T)) { 
  
  print(f)
  load(f)

  exp_cons_all[[f]] <- exp_cons
}


do.call(bind_rows, exp_cons_all) %>% filter(num_count == 1 | is_exp) %>% filter(min_exp == 10) %>% filter(AlleleType == "Ins") %>% filter(sum_count >= 3) %>% mutate(cor = (sum_count == max_count)) %>% group_by(AlleleType, num_count, is_exp, cor, min_prob) %>% summarise(num = sum(n)) %>% print(n = 200)


exp_count_min1 <- do.call(bind_rows, exp_count_all) %>%
  add_column(min_num_exp = ">0")
  
exp_count_min2 <- do.call(bind_rows, exp_count_all) %>%
  filter(num_exp >= 2) %>%
  add_column(min_num_exp = ">1")

var_cols <- wes_palette("Rushmore1")[c(3,1,4)]

p1 <- exp_count_min1 %>%
  rbind(exp_count_min2) %>%
  filter(AlleleType != "Ref") %>%
  filter(AlleleType != "Com") %>%
  group_by(AlleleType, min_exp, min_num_exp) %>%
  summarise(total = sum(n)) %>%
  mutate(AlleleType = recode_factor(AlleleType, "SNV" = "SNV", "Ins" = "Insertion", "Del" = "Deletion")) %>%
  mutate(min_num_exp = factor(min_num_exp, levels=c(">0", ">1"))) %>%
  ggplot(aes(y = log10(total), x = min_exp, color = AlleleType, linetype = min_num_exp, shape = min_num_exp)) +
  geom_line(size = 0.75) +
  geom_point(size = 1.5) +
  annotation_logticks(sides = "l") +
  scale_y_continuous(breaks = seq(2, 5), labels = c(100, 1000, 10000, 100000)) + 
  scale_color_manual(values = var_cols) +
  xlab("Minimum allele expression (TPM)") +
  ylab("Number of expressed alleles") +
  labs(color = "Allele\ntype") +
  labs(linetype = "Number\nof tissues", shape = "Number\nof tissues") + 
  theme_bw() +
  theme(aspect.ratio = 2) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(0.75, "cm")) +
  guides(color = guide_legend(ncol = 1)) +
  guides(linetype = guide_legend(ncol = 1)) +
  guides(shape = guide_legend(ncol = 1)) +
  theme(legend.position = 'bottom') +
  theme(text = element_text(size = 11))



exp_cons_hpa <- do.call(bind_rows, exp_cons_all) %>%
  add_column(hp_filt = "All")

exp_count_hp5 <- do.call(bind_rows, exp_cons_all) %>%
  filter(HomopolymerLength < 6) %>%
  add_column(hp_filt = "<6")

p2 <- exp_cons_hpa %>%
  rbind(exp_count_hp5) %>%
  filter(num_count == 1 | is_exp) %>%
  filter(sum_count >= 2) %>%
  filter(AlleleType != "Ref") %>%
  filter(AlleleType != "Com") %>%
  group_by(AlleleType, min_exp, hp_filt) %>%
  summarise(total = sum(n), cor = sum((sum_count == max_count) * n)) %>%
  mutate(AlleleType = recode_factor(AlleleType, "SNV" = "SNV", "Ins" = "Insertion", "Del" = "Deletion")) %>%
  mutate(hp_filt = factor(hp_filt, levels=c("All", "<6"))) %>%
  ggplot(aes(y = cor / total, x = min_exp, color = AlleleType, linetype = hp_filt, shape = hp_filt)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  ylim(c(0.5, 1)) +
  scale_color_manual(values = var_cols) +
  xlab("Minimum variant expression (TPM)") +
  ylab("Tissue concordance for alleles\nin expressed exons") +
  labs(color = "Allele type") + 
  labs(linetype = "Homo-\npolymer\nlength", shape = "Homo-\npolymer\nlength") + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.5, "cm")) +
  theme(panel.border = element_rect(size = 1)) +
  theme(legend.key.height = unit(1.0, 'cm')) +
  theme(text = element_text(size = 20))

p3 <- exp_cons_hpa %>%
  rbind(exp_count_hp5) %>%
  filter(count >= 2 & is_exp) %>%
  filter(AlleleType != "Ref") %>%
  filter(AlleleType != "Com") %>%
  group_by(AlleleType, min_exp, hp_filt) %>%
  summarise(total = sum(n), cor = sum((sum_count == max_count) * n)) %>%
  mutate(AlleleType = recode_factor(AlleleType, "SNV" = "SNV", "Ins" = "Insertion", "Del" = "Deletion")) %>%
  mutate(hp_filt = factor(hp_filt, levels=c("All", "<6"))) %>%
  ggplot(aes(y = cor / total, x = min_exp, color = AlleleType, linetype = hp_filt, shape = hp_filt)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  ylim(c(0.5, 1)) +
  scale_color_manual(values = var_cols) +
  xlab("Minimum variant expression (TPM)") +
  ylab("Tissue concordance for alleles\nexpressed in >1 tissues") +
  labs(color = "Allele type") + 
  labs(linetype = "Homo-\npolymer\nlength", shape = "Homo-\npolymer\nlength") + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.5, "cm")) +
  theme(panel.border = element_rect(size = 1)) +
  theme(text = element_text(size = 20))


pdf("plots/cons/variant_r2_consistency_5_tissues_count.pdf", height = 5, width = 3.5, pointsize = 12)
print(p1)
dev.off()

pdf("plots/cons/variant_r2_consistency_5_tissues_cons_exon.pdf", height = 5, width = 7, pointsize = 12)
print(p2)
dev.off()

pdf("plots/cons/variant_r2_consistency_5_tissues_cons_var.pdf", height = 5, width = 7, pointsize = 12)
print(p3)
dev.off()

########

calc_zygo_inner <- function(min_exp, min_prob, exp_data) {
  
  print(min_prob)
  print(min_exp)
  
  exp_cons <- exp_data %>% 
    filter(var_exp >= min_exp) %>%
    filter(Probability >= min_prob & TPM > 0) %>%
    group_by(Chrom, Position, Sample) %>%
    mutate(var_num_exp = n()) %>%
    group_by(Chrom, Position, AlleleNum, AlleleType) %>%
    summarise(num_exp = n(), num_exp_hom = sum(var_num_exp == 1), num_exp_het = sum(var_num_exp == 2)) %>%
    mutate(AlleleType = as.factor(AlleleType)) %>%
    add_column(min_prob = min_prob) %>%
    add_column(min_exp = min_exp)
  
  return(exp_cons)
}

parse_and_analyse_exp_vep_data <- function(chrom) {
  
  gc()
  
  print(chrom)  
  
  exp_data <- map_dfr(list.files(path = paste("allele_expression/", chrom, sep = ""), pattern = "allele_exp_rpvg_mpmap.*.txt", full.names = T, recursive = T), parse_exp_data)
  
  print(unique(exp_data$Sample))
  
  exp_cons <- map_dfr(c(1, 5, 10), calc_zygo_inner, 0.8, exp_data)
  
  vep <- read_table(paste("../variant_r1/effect_prediction/1kg_all_exons_", chrom, "_vep_subset.txt.gz", sep = ""), col_names = T, comment = '##', col_types = paste(rep("c", 35), collapse = "")) %>% 
    separate(Location, c("Chrom", "Pos"), ":", extra = "drop", fill = "right") %>%
    separate(Pos, c("Pos", "Pos2"), "-", extra = "drop", fill = "right") %>%
    mutate(Pos = as.double(Pos)) %>%
    mutate(Pos = ifelse(Allele == "-" & MINIMISED != "1", Pos - 1, Pos)) %>%
    rename(Position = Pos) %>%
    mutate(ALLELE_NUM = as.double(ALLELE_NUM)) %>%
    rename(AlleleNum = ALLELE_NUM)
    
  exp_vep <- exp_cons %>%
    filter(AlleleType != "Ref") %>%
    inner_join(vep, by = c("Chrom", "Position", "AlleleNum"))
  
  save(exp_vep, file = paste("rdata/exp_vep_", chrom, ".RData", sep = "", collapse = ""))
}


map_dfr(c(seq(1, 22), "X"), parse_and_analyse_exp_vep_data)

exp_vep_all <- list()

for (f in list.files(pattern = "exp_vep_.*RData", full.names = T, recursive = T)) { 
  
  print(f)
  load(f)
  
  exp_vep_all[[f]] <- exp_vep
}

do.call(bind_rows, exp_vep_all) %>% filter(is.na(Allele))

vep_exp <- do.call(bind_rows, exp_vep_all) %>%
  filter(min_exp == 5) %>%
  filter(IMPACT == "HIGH" | IMPACT == "MODERATE") %>%
  mutate(Consequence = gsub("NMD_transcript_variant", "", Consequence, fixed = T)) %>%
  mutate(Consequence = gsub("splice_region_variant", "", Consequence, fixed = T)) %>%
  mutate(Consequence = gsub("start_retained_variant", "", Consequence, fixed = T)) %>%
  mutate(Consequence = gsub("stop_retained_variant", "", Consequence, fixed = T)) %>%
  mutate(Consequence = gsub("non_coding_transcript_variant", "", Consequence, fixed = T)) %>%
  mutate(Consequence = gsub("synonymous_variant", "", Consequence, fixed = T)) %>%
  filter(Consequence != ",,") %>%
  filter(Consequence != ",") %>%
  filter(Consequence != "")

vep_exp_var <- vep_exp %>%
  group_by(Chrom, Position, REF_ALLELE, Allele, num_exp_hom, num_exp_het) %>%
  summarise(Consequence = paste0(Consequence, collapse = ","), IMPACT = paste0(IMPACT, collapse = ",")) %>%
  rowwise() %>%
  mutate(Consequence = paste(sort(unique(str_split(Consequence, ",")[[1]])), collapse = " / ")) %>%
  mutate(IMPACT = paste(sort(unique(str_split(IMPACT, ",")[[1]])), collapse = " / ")) %>%
  mutate(Consequence = gsub("^ \\/ ", "", Consequence)) %>%
  group_by(Consequence, num_exp_hom, num_exp_het) %>%
  summarise(n = n()) %>%
  add_column(Feature = "Variant")

vep_exp_txp <- vep_exp %>%
  group_by(Chrom, Position, REF_ALLELE, Allele, Feature, num_exp_hom, num_exp_het) %>%
  summarise(Consequence = paste0(Consequence, collapse = ","), IMPACT = paste0(IMPACT, collapse = ",")) %>%
  rowwise() %>%
  mutate(Consequence = paste(sort(unique(str_split(Consequence, ",")[[1]])), collapse = " / ")) %>%
  mutate(IMPACT = paste(sort(unique(str_split(IMPACT, ",")[[1]])), collapse = " / ")) %>%
  mutate(Consequence = gsub("^ \\/ ", "", Consequence)) %>%
  group_by(Consequence, num_exp_hom, num_exp_het) %>%
  summarise(n = n()) %>%
  add_column(Feature = "Transcript")

vep_exp_var$Consequence = recode_factor(vep_exp_var$Consequence, "protein_altering_variant" = "missense_variant")
vep_exp_txp$Consequence = recode_factor(vep_exp_txp$Consequence, "protein_altering_variant" = "missense_variant")

vep_exp_var$Consequence = recode_factor(vep_exp_var$Consequence, 
                                     "missense_variant / splice_acceptor_variant" = "splice_acceptor_variant", 
                                     "missense_variant / splice_donor_variant" = "splice_donor_variant", 
                                     "missense_variant / start_lost" = "start_lost",
                                     "missense_variant / start_lost / stop_lost" = "start_lost / stop_lost",
                                     "missense_variant / stop_gained" = "stop_gained",
                                     "missense_variant / stop_lost" = "stop_lost")
  

vep_exp_var$Consequence = factor(vep_exp_var$Consequence, levels = rev(str_sort(unique(vep_exp_var$Consequence))))
vep_exp_txp$Consequence = factor(vep_exp_txp$Consequence, levels = rev(str_sort(unique(vep_exp_txp$Consequence))))


pdf("plots/vep/variant_r2_vep_var_nonmis_geno.pdf", height = 4, pointsize = 12)
vep_exp_var %>% 
  filter(Consequence != "missense_variant") %>%
  add_column(Genotype = "Inconsistent") %>% 
  mutate(Genotype = ifelse(num_exp_hom > 0 & num_exp_het == 0, "Homozygote", Genotype)) %>% 
  mutate(Genotype = ifelse(num_exp_hom == 0 & num_exp_het > 0, "Heterozygote", Genotype)) %>%
  group_by(Consequence, Genotype) %>%
  summarise(num = sum(n)) %>%
  ggplot(aes(x = Consequence, y = num, fill = Genotype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Rushmore1")[c(3,1,4)]) +
  coord_flip() +
  xlab("") +
  ylab("Number of variants") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

pdf("plots/vep/variant_r2_vep_var_mis_geno.pdf", height = 4, pointsize = 12)
vep_exp_var %>% 
  filter(Consequence == "missense_variant") %>%
  add_column(Genotype = "Inconsistent") %>% 
  mutate(Genotype = ifelse(num_exp_hom > 0 & num_exp_het == 0, "Homozygote", Genotype)) %>% 
  mutate(Genotype = ifelse(num_exp_hom == 0 & num_exp_het > 0, "Heterozygote", Genotype)) %>%
  group_by(Consequence, Genotype) %>%
  summarise(num = sum(n)) %>%
  ggplot(aes(x = Consequence, y = num, fill = Genotype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Rushmore1")[c(3,1,4)]) +
  coord_flip() +
  xlab("") +
  ylab("Number of variants") +
  theme_bw() +
  theme(aspect.ratio = 1 / (length(unique(vep_exp_var$Consequence)) - 1)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

vep_exp_var_hom_counts <- vep_exp_var %>% 
  filter(num_exp_hom > 0 & num_exp_het == 0) %>%
  group_by(Consequence, num_exp_hom) %>%
  summarise(num = sum(n))

vep_exp_var_hom_counts$num_exp_hom = factor(vep_exp_var_hom_counts$num_exp_hom, levels = rev(sort(unique(vep_exp_var_hom_counts$num_exp_hom))))

pdf("plots/vep/variant_r2_vep_var_nonmis_tissues.pdf", height = 3, width = 6, pointsize = 12)
vep_exp_var_hom_counts %>%
  filter(Consequence != "missense_variant") %>%
  ggplot(aes(x = Consequence, y = num, fill = num_exp_hom)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  coord_flip() +
  xlab("") +
  ylab("Number of homozygote variants") +
  labs(fill = "Number\nof tissues") + 
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

pdf("plots/vep/variant_r2_vep_var_mis_tissues.pdf", height = 3, width = 6, pointsize = 12)
vep_exp_var_hom_counts %>% 
  filter(Consequence == "missense_variant") %>%
  ggplot(aes(x = Consequence, y = num, fill = num_exp_hom)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  coord_flip() +
  xlab("") +
  ylab("Number of homozygote variants") +
  labs(fill = "Number\nof tissues") + 
  theme_bw() +
  theme(aspect.ratio = 1 / (length(unique(vep_exp_var_hom_counts$Consequence)))) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

########
