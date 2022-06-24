
# plot_mapping_memory.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r2/")

########

convertMemoryLine <- function(memory_line) {

  return(as.double(strsplit(memory_line, " ")[[1]][6]) / 10^6)
}

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  base_split <- strsplit(basename(filename), "-")[[1]]
  
  mem <- convertMemoryLine(grep('Maximum', readLines(filename), value = T)[1])
  
  data <- data_frame(Mem = mem)
  data <- data %>%
    add_column(Method = base_split[4]) %>%
    add_column(Reads = base_split[7]) %>%
    add_column(Reference = base_split[8])
  
  if (grepl("-f-", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "fast", sep = "_")) %>%
      mutate(Reads = base_split[8]) %>%
      mutate(Reference = base_split[9])
  }
  
  if (grepl("-gs-", basename(filename))) {
    
    data <- data %>%
      mutate(Reference = paste(Reference, "gs", sep = "_"))
  
  } else if (grepl("-gt10-", basename(filename))) {
  
    data <- data %>%
      mutate(Reference = paste(Reference, "gt10", sep = "_"))
  }
  
  return(data)    
}

memory_data_hisat2 <- map_dfr(list.files(path = "../mapping_r1/methods", pattern="*-hisat2-real-.*log.txt", full.names = T, recursive = T), parse_file) 
memory_data_star <- map_dfr(list.files(path = "../mapping_r1/methods", pattern=".*star-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_map_fast <- map_dfr(list.files(path = "./methods", pattern=".*-map-f-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_mpmap <- map_dfr(list.files(path = "./methods", pattern=".*-mpmap-real-.*log.txt", full.names = T, recursive = T), parse_file)

memory_data <- bind_rows(memory_data_hisat2, memory_data_star, memory_data_map_fast, memory_data_mpmap)

########

memory_data$Method = recode_factor(memory_data$Method, 
                                    "hisat2" = "HISAT2", 
                                    "star" = "STAR", 
                                    "map_fast" = "vg map", 
                                    "mpmap" = "vg mpmap")

memory_data$Reads = recode_factor(memory_data$Reads, 
                                         "aed1" = "ENCSR000AED_rep1", 
                                         "t2t1" = "CHM13_rep1", 
                                         "470" = "SRR1153470")

memory_data <- memory_data %>%
  filter(Reference != "na") %>%
  filter(Reference != "eurnceu")

for (reads in unique(memory_data$Reads)) {
  
  memory_data_reads <- memory_data %>%
    filter(Reads == reads)
  
  if (reads == "ENCSR000AED_rep1" | reads == "SRR1153470") {
    
    memory_data_reads <- memory_data_reads %>%
      filter(Reference != "all") %>%
      add_row(Mem = 0, Method = "STAR", Reference = "nceu")
  }
  
  if (reads == "CHM13_rep1") {
    
    memory_data_reads <- memory_data_reads %>%
      add_row(Mem = 0, Method = "STAR", Reference = "all")
  }
  
  memory_data_reads$Reference = recode_factor(memory_data_reads$Reference, 
                                                 "gc100" = "Spliced\nreference",
                                                 "nceu" = "Spliced\npangenome\ngraph",
                                                 "all" = "Spliced\npangenome\ngraph")
  
  memory_data_reads$FacetCol <- "Real reads"
  memory_data_reads$FacetRow <- ""
  
  plotMappingMemoryBenchmark(memory_data_reads, wes_cols, paste("plots/real_memory/real_r2_mapping_memory_", reads, sep = ""))
}

########
