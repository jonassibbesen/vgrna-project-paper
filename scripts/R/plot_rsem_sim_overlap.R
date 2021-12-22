
# plot_rsem_sim_overlap.R

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
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

overlap_data_o1 <- map_dfr(list.files(pattern=".*_ol_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o1$Threshold <- "Overlap >= 1%"

overlap_data_o1 <- overlap_data_o1 %>%
  mutate(Correct = Overlap > 0.01)

overlap_data_o50 <- map_dfr(list.files(pattern=".*_ol_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o50$Threshold <- "Overlap >= 50%"

overlap_data_o50 <- overlap_data_o50 %>%
  mutate(Correct = Overlap > 0.5)

overlap_data_o90 <- map_dfr(list.files(pattern=".*_ol_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o90$Threshold <- "Overlap >= 90%"

overlap_data_o90 <- overlap_data_o90 %>%
  mutate(Correct = Overlap > 0.90)

# overlap_data_o90 %>%
#   filter(Reads == "SRR1153470_uni") %>%
#   filter(Method == "mpmap") %>%
#   filter(Graph == "1kg_NA12878_gencode100") %>%
#   filter(MapQ == 60) %>%
#   ggplot(aes(x = SoftClipLength, color = Correct)) +
#   geom_histogram(fill="white", position="dodge")

overlap_data <- bind_rows(overlap_data_o90, overlap_data_o50, overlap_data_o1) %>%
  filter(Reads == "SRR1153470_uni") 
  
overlap_data$Threshold <- factor(overlap_data$Threshold, levels = c("Overlap >= 90%", "Overlap >= 50%", "Overlap >= 1%"))
overlap_data$Method = recode_factor(overlap_data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
overlap_data$Reads = recode_factor(overlap_data$Reads, "SRR1153470_uni" = "Training set", "ENCSR000AED_rep1_uni" = "Test set")

overlap_data_all <- overlap_data %>%
  filter(Method != "map_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af01_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes")

overlap_data_all$Graph = recode_factor(overlap_data_all$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")
plotDistanceBenchmark(overlap_data_all, wes_cols, "rsem_sim_benchmark_overlap")

overlap_data_sj <- overlap_data %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_NA12878_exons_gencode100") %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af01_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") 

overlap_data_sj$Graph = recode_factor(overlap_data_sj$Graph, "gencode100" = "All transcripts", "1kg_nonCEU_af001_gencode100" = "All transcripts", "gencode85" = "85% transcripts", "1kg_nonCEU_af001_gencode85" = "85% transcripts")
plotDistanceBenchmark(overlap_data_sj, wes_cols, "rsem_sim_benchmark_overlap_sj")


overlap_data_paths <- overlap_data %>%
  filter(Reads == "Training set") %>%
  filter(Graph == "1kg_nonCEU_af001_gencode100" | Method == "STAR") 

overlap_data_paths[overlap_data_paths$Method == "vg map",]$Graph <- "With transcript paths"
overlap_data_paths[overlap_data_paths$Method == "map_nopaths",]$Graph <- "Without transcript paths"
overlap_data_paths[overlap_data_paths$Method == "map_nopaths",]$Method <- "vg map"

overlap_data_paths[overlap_data_paths$Method == "vg mpmap",]$Graph <- "With transcript paths"
overlap_data_paths[overlap_data_paths$Method == "mpmap_nopaths",]$Graph <- "Without transcript paths"
overlap_data_paths[overlap_data_paths$Method == "mpmap_nopaths",]$Method <- "vg mpmap"

overlap_data_paths[overlap_data_paths$Method == "STAR",]$Graph <- "Without transcript paths"
overlap_data_paths[overlap_data_paths$Method == "HISAT2",]$Graph <- "Without transcript paths"

overlap_data_paths$Graph <- factor(overlap_data_paths$Graph, levels = c("With transcript paths", "Without transcript paths"))
plotDistanceBenchmark(overlap_data_paths, wes_cols, "rsem_sim_benchmark_overlap_paths")
