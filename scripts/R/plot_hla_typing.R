
# plot_hla_typing.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/hla_r2/")

########

inference <- "_r2"
samples <- c("NA06994", "NA07037", "NA07357", "NA11829", "NA11893", "NA12006", "NA12043", "NA12234", "NA12272", "NA12275")

#samples <- c("NA07051", "NA11832", "NA11840", "NA11930", "NA12287")

parse_rpvg <- function(filename) {
  
  print(filename)
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Method = dir_split[3]) %>%
    add_column(Sample = dir_split[4])

  return(data)
}

info <- read_table(gzfile("../hla_r1/graphs/1kg_nonCEU_af001_imgt_hla_p10k_noB258_noN_main_gencode100/6/1kg_nonCEU_af001_imgt_hla_p10k_noB258_noN_main_gencode100_6.txt.gz"), col_names = T) %>%
  mutate(Haplotypes = gsub("_thread_hla_", "", Haplotypes, fixed = T)) %>%
  mutate(Haplotypes = gsub("_0_0_0", ":", Haplotypes, fixed = T)) %>%
  mutate(Haplotypes = gsub("*", ":", Haplotypes, fixed = T)) %>%
  select(-Length, -Transcript, -Reference) %>%
  separate_rows(Haplotypes, sep = ",")

info$Haplotypes <- gsub("^([A-Z,0-9]*:[0-9]*:[0-9]*:[0-9]*):.*", "\\1:", info$Haplotypes)

group <- read_table("../hla_r1/groups/hla_nom_g.txt", comment = "#", col_names = F)
group <- group %>%
  add_column(g_id = seq(1, nrow(group))) %>%
  separate(X1, c("Gene", "Allele"), sep = "\\*") %>%
  mutate(Allele = gsub("^;", "", Allele)) %>%
  mutate(Allele = gsub("/", ";", Allele, fixed = T)) %>%
  mutate(Allele = gsub("[A-Z]$", ":", Allele)) %>%
  mutate(Allele = gsub(";", ":,;:", Allele, fixed = T)) %>%
  mutate(Allele = gsub(",;:$", "", Allele)) %>%
  rowwise() %>%
  mutate(Allele = gsub("^", ":", Allele)) %>%
  mutate(Allele = gsub("^", Gene, Allele)) %>%
  mutate(Allele = gsub(";", Gene, Allele, fixed = T)) %>%
  select(-Gene) %>%
  separate_rows(Allele, sep = ",")

group$Allele <- gsub("^([A-Z,0-9]*:[0-9]*:[0-9]*:[0-9]*):.*", "\\1:", group$Allele)

group2 <- group
group2$Allele <- gsub("^([A-Z,0-9]*:[0-9]*:[0-9]*):.*", "\\1:", group2$Allele)

group <- group %>%
  rbind(group2) %>%
  arrange(Allele, g_id) %>%
  distinct()

info <- info %>%
  left_join(group, by = c("Haplotypes" = "Allele"))

info <- info %>%
  mutate(g_id = ifelse(Haplotypes == "C:04:61:", max(g_id + 1, na.rm = T), g_id)) 

info %>% filter(Haplotypes == "C:04:61:")

info <- info %>%
  arrange(Name, Haplotypes, g_id) %>%
  distinct()

main_hla_genes <- read_table("../hla_r1/genes/gencode.v29.primary_assembly.annotation_renamed_full_gene_transcripts.txt", col_names = F) %>%
  filter(X4 == "HLA-A" | X4 == "HLA-B" | X4 == "HLA-C" | X4 == "HLA-DRB1" | X4 == "HLA-DQB1" | X4 == "HLA-DQB1-AS1")

eval <- function(truth, samples) {

  results <- list()
  
  for (sam in samples) {
    
    truth_sam <- truth %>%
      filter(id == sam) %>%
      gather() %>%
      filter(key != "id") %>%
      separate_rows(value, sep = "/")
    
    truth_sam <- truth_sam %>%
      add_column(is_alt = grepl(".1", truth_sam$key, fixed = T)) %>% 
      mutate(key = gsub(".1", "", key, fixed = T)) %>%
      mutate(key = gsub(".2", "", key, fixed = T)) %>%
      mutate(key = gsub("HLA.", "", key, fixed = T)) %>%
      mutate(key = gsub("*", "", key, fixed = T)) %>%
      mutate(value = gsub("*", "", value, fixed = T)) %>%
      unite(name, key:value, sep = ":")
    
    truth_sam$name <- paste(truth_sam$name, ":", sep = "")
    truth_sam$full <- truth_sam$name
    
    truth_sam$d2 <- truth_sam$name
    truth_sam$d2 <- gsub("^([A-Z,0-9]*:[0-9]*:[0-9]*):.*", "\\1:", truth_sam$d2)
    
    truth_sam$d1 <- truth_sam$name
    truth_sam$d1 <- gsub("^([A-Z,0-9]*:[0-9]*):.*", "\\1:", truth_sam$d1)
    
    truth_sam <- truth_sam %>% 
      left_join(group, by = c("name" = "Allele"))
    
    info_correct <- info %>% 
      mutate(correct_full_1 = grepl(paste(unique(truth_sam[!truth_sam$is_alt,]$full), collapse="|"), Haplotypes)) %>% 
      mutate(correct_full_2 = grepl(paste(unique(truth_sam[truth_sam$is_alt,]$full), collapse="|"), Haplotypes)) %>% 
      mutate(correct_d2_1 = grepl(paste(unique(truth_sam[!truth_sam$is_alt,]$d2), collapse="|"), Haplotypes)) %>% 
      mutate(correct_d2_2 = grepl(paste(unique(truth_sam[truth_sam$is_alt,]$d2), collapse="|"), Haplotypes)) %>% 
      mutate(correct_d1_1 = grepl(paste(unique(truth_sam[!truth_sam$is_alt,]$d1), collapse="|"), Haplotypes)) %>% 
      mutate(correct_d1_2 = grepl(paste(unique(truth_sam[truth_sam$is_alt,]$d1), collapse="|"), Haplotypes)) %>%
      mutate(correct_g_1 = g_id %in% truth_sam[!truth_sam$is_alt,]$g_id) %>% 
      mutate(correct_g_2 = g_id %in% truth_sam[truth_sam$is_alt,]$g_id) %>% 
      separate(Haplotypes, c("Gene"), sep = ":", extra = "drop") %>%
      group_by(Name, Gene) %>%
      summarise(correct_full_1 = any(correct_full_1), correct_full_2 = any(correct_full_2), correct_d2_1 = any(correct_d2_1), correct_d2_2 = any(correct_d2_2), correct_d1_1 = any(correct_d1_1), correct_d1_2 = any(correct_d1_2), correct_g_1 = any(correct_g_1), correct_g_2 = any(correct_g_2))
    
    exp_data <- map_dfr(list.files(path = paste("rpvg/geuvadis/inference", inference, sep = ""), pattern = paste(".*", sam, ".*1kg_nonCEU_af001_imgt_hla_p10k_noB258_noN_main_gencode100.txt.gz", sep = ""), full.names = T, recursive = T), parse_rpvg)
    
    exp_data_filt <- exp_data %>%
      filter(HaplotypeProbability >= 0.8)
    
    exp_data_filt <- exp_data_filt %>% 
      left_join(info_correct, by = c("Name"))  %>%
      separate(Name, c("transcript", "hap_id"), "_")
    
    exp_data_filt <- exp_data_filt %>%
      filter(transcript %in% main_hla_genes$X2)
    
    print(sam)
    exp_data_filt %>% group_by(ClusterID, Gene) %>% summarise(n = n()) %>% print(n = 100)
    
    results[[sam]] <- exp_data_filt
  }
  
  return(results)
}

truth_old <- read.table("../hla_r1/truth/20140702_hla_diversity.txt", header = T)
truth_old$sbgroup <- NULL

truth_new <- read.table("../hla_r1/truth/20181129_HLA_types_full_1000_Genomes_Project_panel.txt", header = T, sep = "\t") %>%
  rename(id = Sample.ID)
truth_new$Region <- NULL
truth_new$Population <- NULL

results_old <- do.call(rbind, eval(truth_old, samples))
results_new <- do.call(rbind, eval(truth_new, samples))

calc_hap_stats <- function(data) {
  
  data <- data %>% 
    filter(TPM > 0) %>%
    mutate(correct = (correct_1 | correct_2)) 
    
    return(data)
}

plot_hap_stats <- function(stats1, stats2, title) {
  p <- stats1 %>%
    full_join(stats2, by = c("Sample", "Method", "Gene", "transcript", "TPM")) %>%
    replace_na(list(correct.x = F, correct.y = F)) %>%
    mutate(correct = (correct.x | correct.y)) %>%
    group_by(Sample, Method, Gene) %>%
    summarise(n = n(), tpm_correct = sum(correct * TPM), tpm_sum = sum(TPM)) %>%
    ggplot(aes(y = tpm_correct / tpm_sum, x = Gene)) + 
    ggtitle(title) +
    geom_bar(position = "stack", stat = "identity", fill = wes_palette("Chevalier1")[c(1)]) +
    facet_wrap(vars(Sample)) +
    ylab("Fraction of expression") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 12))
  print(p)
}

stats_old_full <- results_old %>%
  mutate(correct_1 = correct_full_1) %>%
  mutate(correct_2 = correct_full_2) %>%
  calc_hap_stats()

stats_new_full <- results_new %>%
  mutate(correct_1 = correct_full_1) %>%
  mutate(correct_2 = correct_full_2) %>%
  calc_hap_stats()

stats_old_d2 <- results_old %>%
  mutate(correct_1 = correct_d2_1) %>%
  mutate(correct_2 = correct_d2_2) %>%
  calc_hap_stats()

stats_new_d2 <- results_new %>%
  mutate(correct_1 = correct_d2_1) %>%
  mutate(correct_2 = correct_d2_2) %>%
  calc_hap_stats()

stats_old_g <- results_old %>%
  mutate(correct_1 = correct_g_1) %>%
  mutate(correct_2 = correct_g_2) %>%
  calc_hap_stats()

stats_new_g <- results_new %>%
  mutate(correct_1 = correct_g_1) %>%
  mutate(correct_2 = correct_g_2) %>%
  calc_hap_stats()

stats_old_d1 <- results_old %>%
  mutate(correct_1 = correct_d1_1) %>%
  mutate(correct_2 = correct_d1_2) %>%
  calc_hap_stats()

stats_new_d1 <- results_new %>%
  mutate(correct_1 = correct_d1_1) %>%
  mutate(correct_2 = correct_d1_2) %>%
  calc_hap_stats()

pdf(paste("plots/geu/hla_r2_hap_stats_samples_geu_inf", inference, ".pdf", sep = ""))

plot_hap_stats(stats_old_full, stats_new_full, "full")
plot_hap_stats(stats_old_d2, stats_new_d2, "2 digit")
plot_hap_stats(stats_old_g, stats_new_g, "G groups")
plot_hap_stats(stats_old_d1, stats_new_d1, "1 digit")

dev.off()


stats_d2 <- stats_old_d2 %>%
  full_join(stats_new_d2, by = c("Sample", "Method", "Gene", "transcript", "TPM")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "2 field")

stats_g <- stats_old_g %>%
  full_join(stats_new_g, by = c("Sample", "Method", "Gene", "transcript", "TPM")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "G group")

stats_d1 <- stats_old_d1 %>%
  full_join(stats_new_d1, by = c("Sample", "Method", "Gene", "transcript", "TPM")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "1 field")

pdf(paste("plots/geu/hla_r2_hap_stats_mean_geu_inf", inference, ".pdf", sep = ""), height = 4, width = 5, pointsize = 12)
stats_d2 %>%
  rbind(stats_g) %>%
  rbind(stats_d1) %>%
  group_by(Sample, Resolution, Gene) %>%
  summarise(frac_tpm = sum(correct * TPM) / sum(TPM)) %>%
  group_by(Resolution, Gene) %>%
  summarise(mean_frac_tpm = mean(frac_tpm)) %>%
  ggplot(aes(y = mean_frac_tpm, x = Gene, fill = Resolution)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = wes_palette("Rushmore1")[c(3,1,4)]) +
  scale_y_continuous(limits = c(0.5, 1), oob = rescale_none) +
  xlab("HLA gene") +
  ylab("Mean fraction TPM on correct haplotypes") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()


calc_dip_stats <- function(data) {
  
  data <- data %>% 
    filter(TPM > 0) %>%
    mutate(correct_both = (correct_1 | correct_2)) %>%
    group_by(Sample, Method, Gene, transcript) %>%
    summarise(n = n(), sum_tpm = sum(TPM), correct_1 = any(correct_1), correct_2 = any(correct_2), correct_all = all(correct_both)) %>%
    mutate(correct = ifelse(n == 2, correct_1 & correct_2 & correct_all, correct_all)) %>%
  
  return(data)
}

plot_dip_stats <- function(stats1, stats2, title) {
  p <- stats1 %>%
    full_join(stats2, by = c("Sample", "Method", "Gene", "transcript")) %>%
    replace_na(list(correct.x = F, correct.y = F)) %>%
    mutate(correct = (correct.x | correct.y)) %>%
    group_by(Sample, Method, Gene) %>%
    summarise(n = n(), Yes = sum(correct), No = sum(!correct)) %>%
    gather("Correct", "value", No, Yes) %>%
    mutate(Correct = factor(Correct, levels=c('No', "Yes"))) %>%
    ggplot(aes(fill = Correct, y = value, x = Gene)) + 
    ggtitle(title) +
    geom_bar(position = "stack", stat = "identity") +
    facet_wrap(vars(Sample)) +
    scale_fill_manual(values = wes_palette("Chevalier1")[c(2,1)]) +
    ylab("Number of transcripts") +
    labs(fill = "Correct\ndiplotype") + 
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 12))
  print(p)
}



stats_old_full <- results_old %>%
  mutate(correct_1 = correct_full_1) %>%
  mutate(correct_2 = correct_full_2) %>%
  calc_dip_stats()

stats_new_full <- results_new %>%
  mutate(correct_1 = correct_full_1) %>%
  mutate(correct_2 = correct_full_2) %>%
  calc_dip_stats()

stats_old_d2 <- results_old %>%
  mutate(correct_1 = correct_d2_1) %>%
  mutate(correct_2 = correct_d2_2) %>%
  calc_dip_stats()

stats_new_d2 <- results_new %>%
  mutate(correct_1 = correct_d2_1) %>%
  mutate(correct_2 = correct_d2_2) %>%
  calc_dip_stats()

stats_old_g <- results_old %>%
  mutate(correct_1 = correct_g_1) %>%
  mutate(correct_2 = correct_g_2) %>%
  calc_dip_stats()

stats_new_g <- results_new %>%
  mutate(correct_1 = correct_g_1) %>%
  mutate(correct_2 = correct_g_2) %>%
  calc_dip_stats()

stats_old_d1 <- results_old %>%
  mutate(correct_1 = correct_d1_1) %>%
  mutate(correct_2 = correct_d1_2) %>%
  calc_dip_stats()

stats_new_d1 <- results_new %>%
  mutate(correct_1 = correct_d1_1) %>%
  mutate(correct_2 = correct_d1_2) %>%
  calc_dip_stats()

pdf(paste("plots/geu/hla_r2_dip_stats_samples_geu_inf", inference, ".pdf", sep = ""))

plot_dip_stats(stats_old_full, stats_new_full, "full")
plot_dip_stats(stats_old_d2, stats_new_d2, "2 field")
plot_dip_stats(stats_old_g, stats_new_g, "G groups")
plot_dip_stats(stats_old_d1, stats_new_d1, "1 field")

dev.off()



pdf(paste("plots/geu/hla_r2_dip_stats_num_geu_inf", inference, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
results_old %>%
  filter(TPM > 0) %>%
  group_by(transcript, Sample, Gene) %>%
  summarise(n = n()) %>%
  group_by(Sample, Gene) %>%
  summarise(n2 = n()) %>%
  ggplot(aes(x = Sample, y = n2, fill = Gene)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  coord_flip() +
  xlab("") +
  ylab("Number of expressed transcripts") +
  labs(fill = "HLA gene") + 
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()


stats_d2 <- stats_old_d2 %>%
  full_join(stats_new_d2, by = c("Sample", "Method", "Gene", "transcript")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "2 field")

stats_g <- stats_old_g %>%
  full_join(stats_new_g, by = c("Sample", "Method", "Gene", "transcript")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "G group")

stats_d1 <- stats_old_d1 %>%
  full_join(stats_new_d1, by = c("Sample", "Method", "Gene", "transcript")) %>%
  replace_na(list(correct.x = F, correct.y = F)) %>%
  mutate(correct = (correct.x | correct.y)) %>%
  add_column(Resolution = "1 field")

pdf(paste("plots/geu/hla_r2_dip_stats_mean_geu_inf", inference, ".pdf", sep = ""), height = 4, width = 5, pointsize = 12)
stats_d2 %>%
  rbind(stats_g) %>%
  rbind(stats_d1) %>%
  group_by(Sample, Resolution, Gene) %>%
  summarise(frac_cor = sum(correct) / n()) %>%
  group_by(Resolution, Gene) %>%
  summarise(mean_frac_cor = mean(frac_cor)) %>%
  ggplot(aes(y = mean_frac_cor, x = Gene, fill = Resolution)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = wes_palette("Rushmore1")[c(3,1,4)]) +
  scale_y_continuous(limits = c(0.5, 1), oob = rescale_none) +
  xlab("HLA gene") +
  ylab("Mean fraction of correct diplotype transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

########

exp_data <- map_dfr(list.files(path = paste("rpvg/hgsvc/inference", inference, sep = ""), pattern = ".*1kg_all_af001_imgt_hla_p10k_noB258_noN_a100_gencode100_unidi.txt.gz", full.names = T, recursive = T), parse_rpvg)

genes <- read_table("../hla_r1/genes/gencode.v29.primary_assembly.annotation_renamed_full_gene_transcripts.txt", col_names = F)

exp_data_filt <- exp_data %>%
  filter(HaplotypeProbability >= 0.8) %>%
  separate(Name, c("transcript", "hap_id"), "_") %>% 
  arrange(transcript) %>%
  left_join(genes, by = c("transcript" = "X2")) %>%
  mutate(X4 = gsub("-AS1$", "", X4))
  
exp_data_filt %>% group_by(ClusterID, Sample, X4) %>% filter(grepl("HLA", X4)) %>% summarise(n = n()) %>% arrange(Sample) %>% print(n=200)

cons_trio_hg0051 <- exp_data_filt %>% 
  filter(Sample == "HG00512" | Sample == "HG00513" | Sample == "HG00514") %>%
  group_by(transcript, Sample) %>%
  add_column(n = 1) %>%
  mutate(cumsum_n = cumsum(n), exp_n = sum(TPM > 0), var_exp = sum(TPM)) %>%
  filter(!(Sample == "HG00514" & exp_n == 0)) %>%
  group_by(transcript) %>%
  mutate(min_var_exp = min(var_exp), sum_var_exp = sum(var_exp)) %>%
  select(-Length, -EffectiveLength, -ReadCount, -TPM, -HaplotypeProbability, -exp_n) %>%
  pivot_wider(names_from = cumsum_n, values_from = hap_id) %>%
  mutate(exp_samples = n()) %>%
  rename(hap1 = `1`, hap2 = `2`) %>%
  mutate(hap2 = ifelse(is.na(hap2), hap1, hap2)) %>%
  select(-ClusterID, -var_exp) %>%
  pivot_wider(names_from = Sample, values_from = c(hap1, hap2)) %>%
  filter(!is.na(hap1_HG00514)) %>%
  add_column(cons = "False") %>%
  mutate(cons = ifelse((hap1_HG00514 %in% c(hap1_HG00512, hap2_HG00512)) & ((hap2_HG00514 %in% c(hap1_HG00513, hap2_HG00513))), "True", cons)) %>%
  mutate(cons = ifelse((hap2_HG00514 %in% c(hap1_HG00512, hap2_HG00512)) & ((hap1_HG00514 %in% c(hap1_HG00513, hap2_HG00513))), "True", cons)) %>%
  mutate(cons = ifelse((NA %in% c(hap1_HG00512, hap2_HG00512)) | ((NA %in% c(hap1_HG00513, hap2_HG00513))), "Unknown", cons)) %>%
  mutate(cons = ifelse(min_var_exp == 0, "Unknown", cons)) %>%
  add_column(child = "HG00514")

cons_trio_hg0073 <- exp_data_filt %>% 
  filter(Sample == "HG00731" | Sample == "HG00732" | Sample == "HG00733") %>%
  group_by(transcript, Sample) %>%
  add_column(n = 1) %>%
  mutate(cumsum_n = cumsum(n), exp_n = sum(TPM > 0), var_exp = sum(TPM)) %>%
  filter(!(Sample == "HG00733" & exp_n == 0)) %>%
  group_by(transcript) %>%
  mutate(min_var_exp = min(var_exp), sum_var_exp = sum(var_exp)) %>%
  select(-Length, -EffectiveLength, -ReadCount, -TPM, -HaplotypeProbability, -exp_n) %>%
  pivot_wider(names_from = cumsum_n, values_from = hap_id) %>%
  mutate(exp_samples = n()) %>%
  rename(hap1 = `1`, hap2 = `2`) %>%
  mutate(hap2 = ifelse(is.na(hap2), hap1, hap2)) %>%
  select(-ClusterID, -var_exp) %>%
  pivot_wider(names_from = Sample, values_from = c(hap1, hap2)) %>%
  filter(!is.na(hap1_HG00733)) %>%
  add_column(cons = "False") %>%
  mutate(cons = ifelse((hap1_HG00733 %in% c(hap1_HG00731, hap2_HG00731)) & ((hap2_HG00733 %in% c(hap1_HG00732, hap2_HG00732))), "True", cons)) %>%
  mutate(cons = ifelse((hap2_HG00733 %in% c(hap1_HG00731, hap2_HG00731)) & ((hap1_HG00733 %in% c(hap1_HG00732, hap2_HG00732))), "True", cons)) %>%
  mutate(cons = ifelse((NA %in% c(hap1_HG00731, hap2_HG00731)) | ((NA %in% c(hap1_HG00732, hap2_HG00732))), "Unknown", cons)) %>%
  mutate(cons = ifelse(min_var_exp == 0, "Unknown", cons)) %>%
  add_column(child = "HG00733")

cons_trio_na1923 <- exp_data_filt %>% 
  filter(Sample == "NA19238" | Sample == "NA19239" | Sample == "NA19240") %>%
  group_by(transcript, Sample) %>%
  add_column(n = 1) %>%
  mutate(cumsum_n = cumsum(n), exp_n = sum(TPM > 0), var_exp = sum(TPM)) %>%
  filter(!(Sample == "NA19240" & exp_n == 0)) %>%
  group_by(transcript) %>%
  mutate(min_var_exp = min(var_exp), sum_var_exp = sum(var_exp)) %>%
  select(-Length, -EffectiveLength, -ReadCount, -TPM, -HaplotypeProbability, -exp_n) %>%
  pivot_wider(names_from = cumsum_n, values_from = hap_id) %>%
  mutate(exp_samples = n()) %>%
  rename(hap1 = `1`, hap2 = `2`) %>%
  mutate(hap2 = ifelse(is.na(hap2), hap1, hap2)) %>%
  select(-ClusterID, -var_exp) %>%
  pivot_wider(names_from = Sample, values_from = c(hap1, hap2)) %>%
  filter(!is.na(hap1_NA19240)) %>%
  add_column(cons = "False") %>%
  mutate(cons = ifelse((hap1_NA19240 %in% c(hap1_NA19238, hap2_NA19238)) & ((hap2_NA19240 %in% c(hap1_NA19239, hap2_NA19239))), "True", cons)) %>%
  mutate(cons = ifelse((hap2_NA19240 %in% c(hap1_NA19238, hap2_NA19238)) & ((hap1_NA19240 %in% c(hap1_NA19239, hap2_NA19239))), "True", cons)) %>%
  mutate(cons = ifelse((NA %in% c(hap1_NA19238, hap2_NA19238)) | ((NA %in% c(hap1_NA19239, hap2_NA19239))), "Unknown", cons)) %>%
  mutate(cons = ifelse(min_var_exp == 0, "Unknown", cons)) %>%
  add_column(child = "NA19240")


cons_trio <- cons_trio_hg0051 %>%
  rbind(cons_trio_hg0073, cons_trio_na1923) %>%
  filter(grepl("HLA", X4)) %>%
  add_column(dummy = "") %>%
  filter(sum_var_exp > 0) %>%
  mutate(X4 = gsub("^HLA-", "", X4))

cons_trio$X4 = factor(cons_trio$X4, levels = rev(sort(unique(cons_trio$X4))))
cons_trio$cons = factor(cons_trio$cons, levels = c("True", "False", "Unknown"))


pdf("plots/hgsvc/hla_r2_hgsvc_trio_num.pdf", height = 5, width = 5, pointsize = 12)
cons_trio %>%
  group_by(X4, cons, child, dummy) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = X4, y = n, fill = cons)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,1,5)]) +
  coord_flip() +
  facet_grid(child~dummy) +
  scale_y_continuous(breaks = seq(0, 10, 2)) + 
  xlab("") +
  ylab("Number of expressed transcripts") +
  labs(fill = "Trio\nconcordant") + 
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.2, "cm")) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(text = element_text(size = 12))
dev.off()  
  
  
pdf("plots/hgsvc/hla_r2_hgsvc_trio_frac.pdf", height = 5, width = 3.5, pointsize = 12)
cons_trio %>%
  group_by(X4, child, dummy) %>%
  mutate(sum_sum_var_exp = sum(sum_var_exp)) %>%
  group_by(X4, child, cons, dummy) %>%
  summarise(n = sum(sum_var_exp / sum_sum_var_exp)) %>%
  ggplot(aes(x = X4, y = n, fill = cons)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,1,5)]) +
  coord_flip() +
  facet_grid(child~dummy) +
  xlab("") +
  ylab("TPM fraction") +
  labs(fill = "Trio\nconcordant") + 
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.2, "cm")) +
  theme(axis.text.y = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

########
