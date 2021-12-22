
# plot_mapping_compute.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


compute_data <- tibble(Time = numeric(), Method = character(), Reads = character(), Graph = character())

compute_data_polya <- compute_data %>%
  add_row(Time = 44467, Method = "vg map", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 153600, Method = "vg map", Reads = "ENCSR000AED_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 10819, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 16323, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 3039.6, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 3109.84, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 1291.49, Method = "STAR", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "ENCSR000AED_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 34933, Method = "vg map", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 63799, Method = "vg map", Reads = "CHM13_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 23103, Method = "vg mpmap", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 36381, Method = "vg mpmap", Reads = "CHM13_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 18099.07, Method = "HISAT2", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 18192.48, Method = "HISAT2", Reads = "CHM13_rep1", Graph = "Spliced pan-\ngenome graph") %>%
  add_row(Time = 1948.83, Method = "STAR", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "CHM13_rep1", Graph = "Spliced pan-\ngenome graph")

parse_ovl_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    summarise(num_reads = sum(Count))
  
  return(data)
}

for (reads in unique(compute_data_polya$Reads)) {
  
  compute_data_polya_reads <- compute_data_polya %>%
    filter(Reads == reads)

  overlap_data <- map_dfr(list.files(path = "./methods", pattern = paste(".*real_", reads, "_exon_ovl_gc.txt.gz", sep = ""), full.names = T, recursive = T), parse_ovl_file) 
  
  compute_data_polya_reads <- compute_data_polya_reads %>%
    mutate(Time = Time * 16) %>%
    mutate(Time = max(overlap_data) / Time / 2) 
  
  compute_data_polya_reads$FacetCol <- "Real reads"
  compute_data_polya_reads$FacetRow <- ""
  
  print(compute_data_polya_reads)
  
  plotMappingComputeBenchmark(compute_data_polya_reads, wes_cols, paste("plots/polya_rna/real_mapping_compute_polya_", reads, sep = ""))
}

########
