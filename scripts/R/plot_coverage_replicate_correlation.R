# plot_coverage_replicate_correlation.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_types = "iiciii")
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

coverage_data_rep1 <- map_dfr(list.files(pattern=".*_real_.*ENCSR958UOC_rep1_exon_cov_mb.txt", full.names = T, recursive = T), parse_file)
coverage_data_rep2 <- map_dfr(list.files(pattern=".*_real_.*ENCSR958UOC_rep2_exon_cov_mb.txt", full.names = T, recursive = T), parse_file)

coverage_data_mq_cor_list <- list()

for (i in c(0, 1, seq(2, 60, 2), 255)) { 
  
  print(i)
  
  coverage_data_rep1_mq <- coverage_data_rep1 %>%
    mutate(Count = ifelse(MapQ < i, 0, Count)) %>%
    mutate(ReadCoverage = Count * ReadCoverage) %>%
    mutate(BaseCoverage = Count * BaseCoverage) %>%
    group_by(AllelePosition, ExonSize, Reads, Method, Graph) %>%
    summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))
  
  coverage_data_rep2_mq <- coverage_data_rep2 %>%
    mutate(Count = ifelse(MapQ < i, 0, Count)) %>%
    mutate(ReadCoverage = Count * ReadCoverage) %>%
    mutate(BaseCoverage = Count * BaseCoverage) %>%
    group_by(AllelePosition, ExonSize, Reads, Method, Graph) %>%
    summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))  
  
  coverage_data_mq <- right_join(coverage_data_rep1_mq, coverage_data_rep2_mq, by = c("AllelePosition", "ExonSize", "Method", "Graph")) %>%
    mutate(BaseCoverage.x_norm = BaseCoverage.x / ExonSize) %>%
    mutate(BaseCoverage.y_norm = BaseCoverage.y / ExonSize)
  
  coverage_data_mq_cor_pear <- coverage_data_mq %>%
    group_by(Method, Graph) %>%
    summarise(sens = sum((BaseCoverage.x + BaseCoverage.y) / ExonSize), cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Pearson")
  
  coverage_data_mq_cor_log_pear <- coverage_data_mq %>%
    group_by(Method, Graph) %>%
    summarise(sens = sum((BaseCoverage.x + BaseCoverage.y) / ExonSize), cor = cor(log(BaseCoverage.x_norm + 1), log(BaseCoverage.y_norm + 1), method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "LogPearson")
  
  coverage_data_mq_cor_spea <- coverage_data_mq %>%
    group_by(Method, Graph) %>%
    summarise(sens = sum((BaseCoverage.x + BaseCoverage.y) / ExonSize), cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "spearman")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Spearman")
  
  coverage_data_mq_cor_exp <- coverage_data_mq %>%
    group_by(Method, Graph) %>%
    summarise(sens = sum((BaseCoverage.x + BaseCoverage.y) / ExonSize), cor = mean((BaseCoverage.x > 0) == (BaseCoverage.y > 0))) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Expressed")
  
  coverage_data_mq_cor_list[[as.character(i)]] <- rbind(coverage_data_mq_cor_pear, coverage_data_mq_cor_log_pear, coverage_data_mq_cor_spea, coverage_data_mq_cor_exp)
}

coverage_data_mq_cor_data <- do.call(rbind, coverage_data_mq_cor_list) %>%
  mutate(sens = sens / (6606824 + 6599671)) %>%
  filter(!is.na(cor))

coverage_data_mq_cor_data$Method <- recode_factor(coverage_data_mq_cor_data$Method, 
                                      "bowtie2" = "Bowtie2",
                                      "bowtie2_vs_end" = "Bowtie2 (vs end)",
                                      "bowtie2_vs_local" = "Bowtie2 (vs local)",
                                      "hisat2_nosplice" = "HISAT2", 
                                      "star_nosplice" = "STAR", 
                                      "star_encode" = "STAR (encode)",
                                      "map" = "vg map (def)", 
                                      "map_fast" = "vg map",
                                      "mpmap" = "vg mpmap (multi)", 
                                      "mpmap_nomulti" = "vg mpmap")

coverage_data_mq_cor_data <- coverage_data_mq_cor_data %>%
  filter(Method != "Bowtie2 (vs end)") %>%
  filter(Method != "Bowtie2 (vs local)") %>%
  filter(Method != "STAR") %>% 
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (multi)")

coverage_data_mq_cor_data$Graph = recode_factor(coverage_data_mq_cor_data$Graph, 
                                       "linear" = "Linear",
                                       "gencode100" = "Linear",
                                       "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

wes_cols_mir_main <- wes_cols[c(6, seq(1, 5))]

coverage_data_mq_cor_data <- coverage_data_mq_cor_data %>%
  filter(cor_type != "Spearman") %>%
  filter(cor_type != "Expressed") %>%
  filter(cor_type != "LogPearson")

pdf("plots/micro_rna/real_cov_cor_rep_mir.pdf", height = 5, width = 7, pointsize = 12, useDingbats = F)
coverage_data_mq_cor_data %>%
  ggplot(aes(y = cor, x = sens, color = Method, linetype = Graph, shape = Graph, label = Threshold)) +
  geom_line(size = 1.5) + 
  geom_point(data = subset(coverage_data_mq_cor_data, Threshold == 0 | Threshold == 1 | (Threshold == 42 & grepl("Bowtie2", Method)) | Threshold == 60 | Threshold == 255), size = 1.5, alpha = 1) +
  geom_text_repel(data = subset(coverage_data_mq_cor_data, Threshold == 0 | Threshold == 1 | (Threshold == 42 & grepl("Bowtie2", Method)) | Threshold == 60| Threshold == 255), size = 2.5, fontface = 2) +     
  scale_color_manual(values = wes_cols_mir_main) +
  xlab("Fraction mapped bases overlapping microRNAs") +
  ylab("Replicate microRNA coverage correlation") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=16))
dev.off()