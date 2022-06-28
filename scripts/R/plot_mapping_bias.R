
# plot_mapping_bias.R

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
  
  data <- read_table(filename, col_types = "iiiciciii")
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Reference = dir_split[9])
 
  if (grepl("sim_vg", dir_split[6])) {
    
    data <- data %>%
      add_column(Simulation = "Simulated reads (vg)")
    
  } else if (grepl("sim_rsem", dir_split[6])) {
    
    data <- data %>%
      add_column(Simulation = "Simulated reads (RSEM)")
    
  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

min_mapq = 30

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_allele_cov.txt.gz", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Method != "map_fast_multi10") %>%
  filter(Reference != "gencode80") %>%
  filter(Reference != "1kg_nonCEU_af001_gencode80") 

coverage_data_mq <- coverage_data %>%
  mutate(UpReadCount = ifelse(MapQ < min_mapq, 0, UpReadCount)) %>%
  mutate(DownReadCount = ifelse(MapQ < min_mapq, 0, DownReadCount)) %>%
  group_by(VariantPosition, AlleleId, AlleleType, RelativeAlleleLength, Type, Reads, Method, Reference, Simulation) %>%
  summarise(UpReadCount = sum(UpReadCount), DownReadCount = sum(DownReadCount)) 

coverage_data_mq <- full_join(coverage_data_mq[coverage_data_mq$AlleleId == 1,], coverage_data_mq[coverage_data_mq$AlleleId == 2,], by = c("VariantPosition", "Type", "Reads", "Method", "Reference", "Simulation"))

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

coverage_data_mq_bias <- coverage_data_mq %>%
  filter(Type == "polya_rna")

coverage_data_mq_bias$Method = recode_factor(coverage_data_mq_bias$Method, 
                                             "hisat2" = "HISAT2",
                                             "hisat2_multi10" = "HISAT2",
                                             "star" = "STAR",
                                             "star_multi10" = "STAR",
                                             "map" = "vg map (def)", 
                                             "map_fast" = "vg map",
                                             "map_fast_multi10" = "vg map",
                                             "mpmap" = "vg mpmap",
                                             "mpmap_multi10" = "vg mpmap",
                                             "star_alleleseq" = "Diploid reference (STAR)",
                                             "star_alleleseq_levio" = "Diploid reference (STAR, LevioSAM)",
                                             "star_wasp" = "WASP (STAR)")

coverage_data_mq_bias[coverage_data_mq_bias$Method == "WASP (STAR)",]$Reference <- "1kg_NA12878_gencode100"

coverage_data_mq_bias$FacetCol <- paste(coverage_data_mq_bias$Simulation, ", primary alignments", sep = "")
coverage_data_mq_bias$FacetRow <- coverage_data_mq_bias$Reference

########

coverage_data_mq_bias_debug <- coverage_data_mq_bias

for (reads in unique(coverage_data_mq_bias_debug$Reads)) {

  coverage_data_mq_bias_debug_reads <- coverage_data_mq_bias_debug %>%
    filter(Reads == reads)

  plotMappingBiasBenchmark(coverage_data_mq_bias_debug_reads, wes_cols, paste("plots/sim_bias/debug/vg_sim_r2_mapping_bias_debug_", reads, sep = ""), 20)
  
  coverage_data_mq_bias_debug_reads$FacetCol <- coverage_data_mq_bias_debug_reads$var
  coverage_data_mq_bias_debug_reads$FacetRow <- paste(coverage_data_mq_bias_debug_reads$Simulation, ",\n primary alignments", sep = "")
  
  plotMappingBiasBinomBenchmark(coverage_data_mq_bias_debug_reads, wes_cols, paste("plots/sim_bias/debug/vg_sim_r2_mapping_bias_binom_a001_debug_", reads, sep = ""), 0.01, 20)
}

########

coverage_data_mq_bias_main <- coverage_data_mq_bias %>%
  filter(Reads == "sim_vg_r2_ENCSR000AED_rep1_uni") %>%
  filter(Method != "WASP (STAR)") %>%
  filter(Method != "Diploid reference (STAR)") %>%
  filter(Method != "Diploid reference (STAR, LevioSAM)") %>%
  filter(Method != "vg map (def)") %>%
  filter(Reference != "1kg_NA12878_gencode100") %>%
  filter(Reference != "1kg_NA12878_exons_gencode100")

coverage_data_mq_bias_main$Reference = recode_factor(coverage_data_mq_bias_main$Reference,
                                            "1kg_nonCEU_af001_gencode100" = "Spliced pan-\ngenome graph",
                                            "gencode100" = "Spliced\nreference")

coverage_data_mq_bias_main$FacetRow <- coverage_data_mq_bias_main$Reference

for (reads in unique(coverage_data_mq_bias_main$Reads)) {

  coverage_data_mq_bias_main_reads <- coverage_data_mq_bias_main %>%
    filter(Reads == reads)

  plotMappingBiasBenchmark(coverage_data_mq_bias_main_reads, wes_cols, paste("plots/sim_bias/vg_sim_r2_mapping_bias_main_", reads, sep = ""), 20)
}

########

coverage_data_mq_bias_binom <- coverage_data_mq_bias %>%
  filter(Reads == "sim_vg_r2_ENCSR000AED_rep1_uni") %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "Diploid reference (STAR)") %>%
  filter(Method != "Diploid reference (STAR, LevioSAM)")

coverage_data_mq_bias_binom$Reference = recode_factor(coverage_data_mq_bias_binom$Reference,
                                                 "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                                 "1kg_NA12878_gencode100" = "Spliced personal graph/reference",
                                                 "1kg_NA12878_exons_gencode100" = "Spliced personal graph/reference",
                                                 "gencode100" = "Spliced reference")

coverage_data_mq_bias_binom$FacetCol <- coverage_data_mq_bias_binom$var
coverage_data_mq_bias_binom$FacetRow <- paste(coverage_data_mq_bias_binom$Simulation, ",\n primary alignments", sep = "")

for (reads in unique(coverage_data_mq_bias_binom$Reads)) {
  
  coverage_data_mq_bias_binom_reads <- coverage_data_mq_bias_binom %>%
    filter(Reads == reads)
  
  plotMappingBiasBinomBenchmark(coverage_data_mq_bias_binom_reads, wes_cols, paste("plots/sim_bias/vg_sim_r2_mapping_bias_binom_a001_", reads, sep = ""), 0.01, 20)
}

coverage_data_mq_bias_binom_nowasp <- coverage_data_mq_bias_binom %>%
  filter(Method != "WASP (STAR)")

for (reads in unique(coverage_data_mq_bias_binom_nowasp$Reads)) {
  
  coverage_data_mq_bias_binom_nowasp_reads <- coverage_data_mq_bias_binom_nowasp %>%
    filter(Reads == reads)
  
  plotMappingBiasBinomBenchmark(coverage_data_mq_bias_binom_nowasp_reads, wes_cols, paste("plots/sim_bias/vg_sim_r2_mapping_bias_binom_a001_nowasp_", reads, sep = ""), 0.01, 20)
}

########
