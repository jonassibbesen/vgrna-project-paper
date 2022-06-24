
# plot_vg_sim_distance.R

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
  
  data <- read_table(filename, col_names = F)
  colnames(data) <- c("Count", "Distance", "IsMapped", "MapQ", "GroupMapQ", "Method")

  data <- data %>%
    select(-Method) %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Reference = dir_split[9])

  if (grepl("_gamp_dist_", basename(filename))) {

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

########

distance_data_raw_h1 <- map_dfr(list.files(path = "./methods", pattern=".*_dist_h1.txt.gz", full.names = T, recursive = T), parse_file) 
distance_data_raw_h2 <- map_dfr(list.files(path = "./methods", pattern=".*_dist_h2.txt.gz", full.names = T, recursive = T), parse_file) 

distance_data <- rbind(distance_data_raw_h1, distance_data_raw_h2)  %>%
  mutate(Correct = Distance <= distance_threshold) %>%
  mutate(MapQ = GroupMapQ)

distance_data$Method <- recode_factor(distance_data$Method, 
                                      "hisat2_multi10" = "HISAT2",
                                      "star_multi10" = "STAR",
                                      "map_fast_multi10" = "vg map", 
                                      "mpmap_multi10" = "vg mpmap (gam)", 
                                      "mpmap_multi10_gamp" = "vg mpmap")

distance_data$Reference = recode_factor(distance_data$Reference, 
                                   "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                   "gencode100" = "Spliced reference")

distance_data$FacetCol <- paste(distance_data$Simulation, ", multi alignments", sep = "")
distance_data$FacetRow <- ""
distance_data$LineType <- distance_data$Reference

########

distance_data_debug <- distance_data

for (reads in unique(distance_data_debug$Reads)) {
  
  distance_data_debug_reads <- distance_data_debug %>%
    filter(Reads == reads)
  
  plotRocBenchmarkMapQ(distance_data_debug_reads, wes_cols, "Reference", paste("plots/sim_distance/debug/vg_sim_r2_distance_debug_dist", distance_threshold, "_", reads, sep = ""))
}

########

distance_data_main <- distance_data %>%
  filter(Method != "vg mpmap (gam)") 

for (reads in unique(distance_data_main$Reads)) {
  
  distance_data_main_reads <- distance_data_main %>%
    filter(Reads == reads)
  
  plotRocBenchmarkMapQ(distance_data_main_reads, wes_cols, "Reference", paste("plots/sim_distance/vg_sim_r2_distance_main_dist", distance_threshold, "_", reads, sep = ""))
}

########
