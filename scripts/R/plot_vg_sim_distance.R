
# plot_vg_sim_distance.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

# source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = F)
  colnames(data) <- c("Count", "Distance", "IsMapped", "MapQ", "Method")

  data <- data %>%
    select(-Method) %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])

  if (grepl("dist_gamp_", basename(filename))) {

    data <- data %>%
      mutate(Method = paste(Method, "gamp", sep = "_"))
  }
  
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

distance_threshold <- 100

distance_data_raw_h1 <- map_dfr(list.files(path = "./methods", pattern=".*_dist_gam.*h1.txt.gz", full.names = T, recursive = T), parse_file) 
distance_data_raw_h2 <- map_dfr(list.files(path = "./methods", pattern=".*_dist_gam.*h2.txt.gz", full.names = T, recursive = T), parse_file) 

distance_data <- rbind(distance_data_raw_h1, distance_data_raw_h2)  %>%
  mutate(Correct = Distance <= distance_threshold) 


########


distance_data <- distance_data %>%
  filter(Type == "polya_rna")

distance_data$Method <- recode_factor(distance_data$Method, 
                                      "hisat2" = "HISAT2",
                                      "star" = "STAR",
                                      "map" = "vg map (def)", 
                                      "map_fast" = "vg map", 
                                      "mpmap" = "vg mpmap (gam)", 
                                      "mpmap_gamp" = "vg mpmap", 
                                      "mpmap_nosplice" = "vg mpmap (gam)",
                                      "mpmap_nosplice_gamp" = "vg mpmap")

distance_data <- distance_data %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (gam)")

distance_data$FacetCol <- distance_data$Simulation
distance_data$FacetRow <- ""


distance_data_main <- distance_data %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes")

distance_data_main$Graph = recode_factor(distance_data_main$Graph, 
                                               "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                               "1kg_NA12878_gencode100" = "Personal reference graph",
                                               "1kg_NA12878_exons_gencode100" = "Personal reference graph",
                                               "gencode100" = "Spliced reference")

for (reads in unique(distance_data_main$Reads)) {
  
  distance_data_main_reads <- distance_data_main %>%
    filter(Reads == reads)
  
  plotRocBenchmarkMapQ(distance_data_main_reads, wes_cols, paste("plots/sim_distance/vg_sim_distance_main_dist", distance_threshold, "_", reads, sep = ""))
}


distance_data_gene <- distance_data %>%
  filter(Graph != "gencode100") %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Method != "STAR") %>%
  filter(Method != "HISAT2")

distance_data_gene$Graph = recode_factor(distance_data_gene$Graph, 
                                         "1kg_nonCEU_af001_gencode100" = "Whole genome graph",
                                         "1kg_nonCEU_af001_gencode100_genes" = "Exons only graph")

for (reads in unique(distance_data_gene$Reads)) {
  
  distance_data_gene_reads <- distance_data_gene %>%
    filter(Reads == reads)
  
  plotRocBenchmarkMapQ(distance_data_gene_reads, wes_cols, paste("plots/sim_distance/vg_sim_distance_gene_dist", distance_threshold, "_", reads, sep = ""))
}


########
