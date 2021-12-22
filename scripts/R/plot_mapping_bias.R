
# plot_mapping_bias.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

#source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  if (grepl("mpmap", filename)) {
    
    data <- read_table2(filename, col_types = "iiiciciii")
    
  } else {
    
    data <- read_table2(filename, col_types = "iiciciii") %>%
      mutate(AllelicMapQ = MapQ)
  }
  
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

min_mapq = 30

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80") 

coverage_data_amq <- coverage_data %>%
  filter(Method == "mpmap") %>%
  mutate(Method = "mpmap_amq") %>%
  mutate(MapQ = AllelicMapQ)

coverage_data <- rbind(coverage_data, coverage_data_amq)

coverage_data_mq <- coverage_data %>%
  mutate(UpReadCount = ifelse(MapQ < min_mapq, 0, UpReadCount)) %>%
  mutate(DownReadCount = ifelse(MapQ < min_mapq, 0, DownReadCount)) %>%
  group_by(VariantPosition, AlleleId, AlleleType, RelativeAlleleLength, Type, Reads, Method, Graph) %>%
  summarise(UpReadCount = sum(UpReadCount), DownReadCount = sum(DownReadCount)) 

coverage_data_mq <- full_join(coverage_data_mq[coverage_data_mq$AlleleId == 1,], coverage_data_mq[coverage_data_mq$AlleleId == 2,], by = c("VariantPosition", "Type", "Reads", "Method", "Graph"))

coverage_data_mq <- coverage_data_mq %>% 
  filter((AlleleType.x != AlleleType.y) & (AlleleType.x == "REF" | AlleleType.y == "REF")) %>% 
  mutate(ref_up = ifelse(AlleleType.x == "REF", UpReadCount.x, UpReadCount.y)) %>%
  mutate(ref_down = ifelse(AlleleType.x == "REF", DownReadCount.x, DownReadCount.y)) %>%
  mutate(alt_up = ifelse(AlleleType.x != "REF", UpReadCount.x, UpReadCount.y)) %>%
  mutate(alt_down = ifelse(AlleleType.x != "REF", DownReadCount.x, DownReadCount.y)) %>%
  mutate(var = ifelse(AlleleType.x == "REF", AlleleType.y, AlleleType.x)) %>%
  mutate(len = ifelse(AlleleType.x == "REF", RelativeAlleleLength.y, RelativeAlleleLength.x)) %>%
  filter(var != "COM") 

coverage_data_mq$var <- factor(coverage_data_mq$var, levels = c("SNV", "INS", "DEL"))

coverage_data_mq$var = recode_factor(coverage_data_mq$var, 
                                     "SNV" = "SNV", 
                                     "INS" = "Insertion", 
                                     "DEL" = "Deletion")


########


coverage_data_mq_bias <- coverage_data_mq %>%
  filter(Type == "polya_rna")

coverage_data_mq_bias$Method = recode_factor(coverage_data_mq_bias$Method, 
                                        "hisat2" = "HISAT2", 
                                        "star" = "STAR", 
                                        "star_wasp" = "STAR (WASP)",
                                        "map" = "vg map (def)", 
                                        "map_fast" = "vg map",
                                        "mpmap" = "vg mpmap", 
                                        "mpmap_amq" = "vg mpmap (amq)")

coverage_data_mq_bias <- coverage_data_mq_bias %>%
  filter(Method != "vg map (def)")

coverage_data_mq_bias$Graph = recode_factor(coverage_data_mq_bias$Graph,
                                        "1kg_nonCEU_af001_gencode100" = "Spliced\npangeome\ngraph",
                                        "1kg_NA12878_gencode100" = "Personal\nreference\ngraph",
                                        "1kg_NA12878_exons_gencode100" = "Personal\nreference\ngraph",
                                        "gencode100" = "Spliced\nreference")

coverage_data_mq_bias$FacetCol <- "Simulated reads"
coverage_data_mq_bias$FacetRow <- coverage_data_mq_bias$Graph

for (reads in unique(coverage_data_mq_bias$Reads)) {
  
  coverage_data_mq_bias_reads <- coverage_data_mq_bias %>%
    filter(Reads == reads)
  
  plotMappingBiasBenchmark(coverage_data_mq_bias_reads, wes_cols, paste("plots/sim_bias/vg_sim_mapping_bias_main_", reads, sep = ""))
}


min_count <- 10

coverage_data_mq_bias_binom <- coverage_data_mq_bias %>%
  ungroup() %>%
  filter(Reads == "SRR1153470_uni") %>%
  mutate(ref = (ref_up + ref_down) / 2) %>%
  mutate(alt = (alt_up + alt_down) / 2) %>%
  mutate(count = ref + alt) %>%
  rowwise() %>%
  mutate(binom_test = binom.test(x = c(as.integer(ref), as.integer(alt)), alternative = "two.sided")$p.value) %>%
  group_by(Reads, Method, Graph, var, count) %>%
  summarise(n = n(), n_binom = sum(binom_test < 0.01))

coverage_data_mq_bias_binom$FacetCol <- coverage_data_mq_bias_binom$var
coverage_data_mq_bias_binom$FacetRow <- "Simulated reads"

coverage_data_mq_bias_binom$Graph = recode_factor(coverage_data_mq_bias_binom$Graph,
                                            "Spliced\npangeome\ngraph" = "Spliced pangeome graph",
                                            "Personal\nreference\ngraph" = "Personal reference graph",
                                            "Spliced\nreference" = "Spliced reference")

pdf("plots/sim_bias/test.pdf", height = 5, width = 9, pointsize = 10)
coverage_data_mq_bias_binom %>%
  filter(count >= min_count) %>%
  group_by(Reads, Method, Graph, var, FacetCol, FacetRow) %>%
  summarise(n = sum(n), n_binom = sum(n_binom)) %>%
  ggplot(aes(y = n_binom / n, x = n, color = Method,shape = Graph)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = wes_cols) +
  facet_grid(FacetRow ~ FacetCol, scales="free_x") +
  xlab("Number of variants (coverage >= 10)") +
  ylab("Fraction variants with binomial test p-value < 0.01") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1, "cm")) +
  theme(text = element_text(size = 8)) 
dev.off()

