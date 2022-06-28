
# plot_coverage_correlation.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r2/")

########

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]

  data <- read_table(filename) 
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Reference = dir_split[9])

  return(data)
}

pb_coverage <- read_table("../mapping_r1/alignments/ENCSR706ANY/ENCSR706ANY_mq30_exon_cov_bam.txt.gz", col_types = "iiciii")
pb_coverage <- pb_coverage %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_real_.*cov_ENCSR706ANY_mq30.txt.gz", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  mutate(MapQ = ifelse(MapQ > 60, 60, MapQ))

coverage_data <- coverage_data %>%
  filter(Type == "polya_rna")

coverage_data$Method <- recode_factor(coverage_data$Method, 
                                                       "hisat2" = "HISAT2", 
                                                       "star" = "STAR",
                                                       "map_fast" = "vg map",
                                                       "mpmap" = "vg mpmap", 
                                                       "star_alleleseq" = "Diploid reference (STAR)",
                                                       "star_alleleseq_levio" = "Diploid reference (STAR, LevioSAM)")

coverage_data$Reference = recode_factor(coverage_data$Reference, 
                                        "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                        "1kg_NA12878_gencode100" = "Spliced personal graph/reference",
                                        "1kg_NA12878_exons_gencode100" = "Spliced personal graph/reference",
                                        "gencode100" = "Spliced reference")

########

coverage_data_pb_mq_corr_list <- list()

for (i in c(0, 1, seq(5, 60, 5))) { 
  
  print(i)
  
  coverage_data_mq <- coverage_data %>%
    mutate(Count = ifelse(MapQ < i, 0, Count)) %>%
    mutate(ReadCoverage = Count * ReadCoverage) %>%
    mutate(BaseCoverage = Count * BaseCoverage) %>%
    group_by(AllelePosition, ExonSize, Type, Reads, Method, Reference) %>%
    summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))
  
  coverage_data_pb_mq <- right_join(pb_coverage, coverage_data_mq, by = c("AllelePosition", "ExonSize")) %>%
    mutate(BaseCoverage.x_norm = BaseCoverage.x / ExonSize) %>%
    mutate(BaseCoverage.y_norm = BaseCoverage.y / ExonSize)
  
  coverage_data_pb_mq_corr_pear <- coverage_data_pb_mq %>%
    group_by(Type, Reads, Method, Reference) %>%
    summarise(Corr = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) 
    
  coverage_data_pb_mq_corr_list[[as.character(i)]] <- coverage_data_pb_mq_corr_pear
}

coverage_data_pb_mq_corr <- do.call(rbind, coverage_data_pb_mq_corr_list)

coverage_data_pb_mq_corr$FacetCol <- "Real reads, primary alignments"
coverage_data_pb_mq_corr$FacetRow <- ""

######## 

coverage_data_pb_mq_corr_main <- coverage_data_pb_mq_corr %>%
  filter(Reference != "Spliced personal graph/reference")

for (reads in unique(coverage_data_pb_mq_corr_main$Reads)) {

  coverage_data_pb_mq_corr_main_reads <- coverage_data_pb_mq_corr_main %>%
    filter(Reads == reads) %>%
    rename(MapQ = Threshold)

  plotIsoSeqCorrelationBenchmark(coverage_data_pb_mq_corr_main_reads, wes_cols, paste("plots/real_corr/real_r2_cov_corr_main_", reads, sep = ""))
}

########

coverage_data_pb_mq_corr_personal <- coverage_data_pb_mq_corr %>%
  filter(Reads == "ENCSR000AED_rep1") %>%
  filter(Method == "STAR" | ((Method == "vg mpmap" | Method == "Diploid reference (STAR)" | Method == "Diploid reference (STAR, LevioSAM)") & Reference == "Spliced personal graph/reference"))

for (reads in unique(coverage_data_pb_mq_corr_personal$Reads)) {
  
  coverage_data_pb_mq_corr_personal_reads <- coverage_data_pb_mq_corr_personal %>%
    filter(Reads == reads) %>%
    rename(MapQ = Threshold)
  
  plotIsoSeqCorrelationBenchmark(coverage_data_pb_mq_corr_personal_reads, wes_cols[c(2,4,5,6)], paste("plots/real_corr/real_r2_cov_corr_personal_", reads, sep = ""))
}

########

coverage_data_mq30 <- coverage_data %>%
  filter(Reads == "ENCSR000AED_rep1") %>%
  filter(Reference != "Spliced personal graph/reference") %>%
  mutate(Count = ifelse(MapQ < 30, 0, Count)) %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize, Type, Reads, Method, Reference) %>%
  summarise(BaseCoverage = sum(BaseCoverage))

coverage_data_pb_mq30 <- right_join(pb_coverage, coverage_data_mq30, by = c("AllelePosition", "ExonSize"), suffix = c(".est", ".pb")) %>%
  mutate(Coverage.est = BaseCoverage.est / ExonSize) %>%
  mutate(Coverage.pb = BaseCoverage.pb / ExonSize)

coverage_data_pb_mq30$Reference = recode_factor(coverage_data_pb_mq30$Reference,
                                    "Spliced pangenome graph" = "Spliced pan-\ngenome graph",
                                    "Spliced reference" = "Spliced\nreference")

coverage_data_pb_mq30$FacetCol <- coverage_data_pb_mq30$Method
coverage_data_pb_mq30$FacetRow <- coverage_data_pb_mq30$Reference

for (reads in unique(coverage_data_pb_mq30$Reads)) {

  coverage_data_pb_mq30_reads <- coverage_data_pb_mq30 %>%
    filter(Reads == reads)

  plotIsoSeqCoverageBenchmark(coverage_data_pb_mq30_reads, wes_cols, paste("plots/real_corr/real_r2_cov_scatter_main_", reads, sep = ""))
}

########
