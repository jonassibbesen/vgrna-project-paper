
# plot_region_overlap.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  if (grepl("ovl_ENCSR706ANY_mq0", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "ENCSR706ANY")
    
  } else if (grepl("ovl_ENCSR706ANY_mq30", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "ENCSR706ANY (MapQ >= 30)")

  } else if (grepl("ovl_gc", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "GENCODE")

  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

overlap_data_raw <- map_dfr(list.files(path = "./methods", pattern=".*_exon_ovl_.*.txt", full.names = T, recursive = T), parse_file)

overlap_data_o05 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.05) %>%
  add_column(Threshold = "Overlap >= 5%") 

overlap_data_o80 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.80) %>%
  add_column(Threshold = "Overlap >= 80%") 

overlap_data <- bind_rows(overlap_data_o05, overlap_data_o80)


########

#dataset <- "CHM13_rep1"
dataset <- "SRR1153470"

overlap_data_polya <- overlap_data %>%
  filter(Reads == dataset) 

overlap_data_polya$Method <- recode_factor(overlap_data_polya$Method, 
                                           "hisat2" = "HISAT2",
                                           "star" = "STAR",
                                           "map_fast" = "vg map", 
                                           "mpmap" = "vg mpmap")

overlap_data_polya$Graph = recode_factor(overlap_data_polya$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "Spliced variant graph",
                                         "1kg_all_af001_gencode100" = "Spliced variant graph")

overlap_data_polya_main <- overlap_data_polya

overlap_data_polya_main$FacetCol <- recode_factor(overlap_data_polya_main$Threshold, 
                                       "Overlap >= 5%" = "Overlap >= 5%",
                                       "Overlap >= 80%" = "Overlap >= 80%")

overlap_data_polya_main$FacetRow <- recode_factor(overlap_data_polya_main$Regions, 
                                       "ENCSR706ANY" = "ENCSR706ANY", 
                                       "ENCSR706ANY (MapQ >= 30)" = "ENCSR706ANY (MapQ >= 30)", 
                                       "GENCODE" = "GENCODE")

plotOverlapBenchmarkMapQ(overlap_data_polya_main, wes_cols, paste("plots/polya_rna/real_overlap_polya_main_", dataset, sep = ""))

