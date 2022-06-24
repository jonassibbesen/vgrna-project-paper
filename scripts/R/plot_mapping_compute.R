
# plot_mapping_compute.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r2/")

########

compute_data <- tibble(Time = numeric(), Method = character(), Reads = character(), Reference = character())

compute_data <- compute_data %>%
  add_row(Time = 44232, Method = "vg map", Reads = "ENCSR000AED_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 146693, Method = "vg map", Reads = "ENCSR000AED_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 9553, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 14262, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 3040, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 3110, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 1306, Method = "STAR", Reads = "ENCSR000AED_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "ENCSR000AED_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 34865, Method = "vg map", Reads = "CHM13_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 62966, Method = "vg map", Reads = "CHM13_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 15752, Method = "vg mpmap", Reads = "CHM13_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 22188, Method = "vg mpmap", Reads = "CHM13_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 18099, Method = "HISAT2", Reads = "CHM13_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = 18192, Method = "HISAT2", Reads = "CHM13_rep1", Reference = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 1963, Method = "STAR", Reads = "CHM13_rep1", Reference = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "CHM13_rep1", Reference = "Spliced\npangenome\ngraph")

parse_ovl_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table(filename)
  data <- data %>%
    summarise(num_reads = sum(Count))
  
  return(data)
}

compute_data$Reference <- factor(compute_data$Reference, levels = c("Spliced\nreference", "Spliced\npangenome\ngraph"))

for (reads in unique(compute_data$Reads)) {
  
  compute_data_reads <- compute_data %>%
    filter(Reads == reads)

  overlap_data <- map_dfr(list.files(path = "./methods", pattern = paste(".*real_r2_", reads, "_exon_ovl_gc.txt.gz", sep = ""), full.names = T, recursive = T), parse_ovl_file) 
  
  compute_data_reads <- compute_data_reads %>%
    mutate(Time = Time * 16) %>%
    mutate(Time = max(overlap_data) / Time / 2) 
  
  compute_data_reads$FacetCol <- "Real reads"
  compute_data_reads$FacetRow <- ""
  
  plotMappingComputeBenchmark(compute_data_reads, wes_cols, paste("plots/real_compute/real_r2_mapping_compute_polya_", reads, sep = ""))
}

########
