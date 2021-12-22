
# plot_mapping_compute_kube.R

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


convertTimeLine <- function(time_line) {
  
  time <- 0
  time_idx = 0
  
  for (value in rev(strsplit(strsplit(time_line, " ")[[1]][8], ":")[[1]])) {
    
    time = time + as.double(value) * 60^time_idx
    time_idx = time_idx + 1
  }
  
  return(time)
}

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  base_split <- strsplit(basename(filename), "-")[[1]]
  
  time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1])

  if (base_split[4] == "hisat2") {
    
    time <- time + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[2])
  }
  
  data <- data_frame(Time = time)
  data <- data %>%
    add_column(Method = base_split[4]) %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = base_split[6]) %>%
    add_column(Graph = base_split[7])

  if (grepl("-f-", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "fast", sep = "_")) %>%
      mutate(Reads = base_split[7]) %>%
      mutate(Graph = base_split[8])
  }
  
  if (grepl("-gs-", basename(filename))) {
    
    data <- data %>%
      mutate(Graph = paste(Graph, "gs", sep = "_"))
  
    } else if (grepl("-gt10-", basename(filename))) {
    
    data <- data %>%
      mutate(Graph = paste(Graph, "gt10", sep = "_"))
  }
  
  return(data)    
}

compute_data_hisat2 <- map_dfr(list.files(pattern="-hisat2-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_star <- map_dfr(list.files(pattern=".*star-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_map_fast <- map_dfr(list.files(pattern=".*-map-f-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-.*log.txt", full.names = T, recursive = T), parse_file)

compute_data <- bind_rows(compute_data_hisat2, compute_data_star, compute_data_map_fast, compute_data_mpmap)


########


compute_data_polya <- compute_data %>%
  filter(Type == "polya_rna")

compute_data_polya$Method = recode_factor(compute_data_polya$Method, 
                                   "hisat2" = "HISAT2", 
                                   "star" = "STAR", 
                                   "map_fast" = "vg map", 
                                   "mpmap" = "vg mpmap")

compute_data_polya$Reads = recode_factor(compute_data_polya$Reads, 
                                    "aed1" = "ENCSR000AED_rep1", 
                                    "t2t1" = "CHM13_rep1", 
                                    "470" = "SRR1153470")

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
  
  if (reads == "ENCSR000AED_rep1" | reads == "SRR1153470") {
    
    compute_data_polya_reads <- compute_data_polya_reads %>%
      filter(Graph != "all") %>%
      add_row(Time = 0, Method = "STAR", Graph = "nceu")
  }
  
  if (reads == "CHM13_rep1") {
    
    compute_data_polya_reads <- compute_data_polya_reads %>%
      add_row(Time = 0, Method = "STAR", Graph = "all")
  }

  compute_data_polya_reads$Graph = recode_factor(compute_data_polya_reads$Graph, 
                                   "gc100" = "Spliced\nreference",
                                   "nceu" = "Spliced pan-\ngenome graph",
                                   "all" = "Spliced pan-\ngenome graph")

  compute_data_polya_reads$FacetCol <- "Real reads"
  compute_data_polya_reads$FacetRow <- ""
  
  print(compute_data_polya_reads)
  
  plotMappingComputeBenchmark(compute_data_polya_reads, wes_cols, paste("plots/polya_rna/real_mapping_compute_polya_kube_", reads, sep = ""))
}

########
