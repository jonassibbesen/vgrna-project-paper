
# plot_mapping_stats.R

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

mapping_data <- map_dfr(list.files(path = "./methods", pattern=".*_exon_ovl_gc.*.txt", full.names = T, recursive = T), parse_file)

mapping_data <- mapping_data %>%
  filter(Type == "polya_rna")

mapping_data$Method <- recode_factor(mapping_data$Method, 
                                     "hisat2" = "HISAT2",
                                     "star" = "STAR",
                                     "map_fast" = "vg map", 
                                     "mpmap" = "vg mpmap",
                                     "star_alleleseq" = "Diploid reference (STAR)",
                                     "star_alleleseq_levio" = "Diploid reference (STAR, LevioSAM)")


mapping_data$Reference = recode_factor(mapping_data$Reference, 
                                            "gencode100" = "Spliced reference",
                                            "1kg_NA12878_gencode100" = "Spliced personal graph/reference",
                                            "1kg_NA12878_exons_gencode100" = "Spliced personal graph/reference",
                                            "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                            "1kg_all_af001_gencode100" = "Spliced pangenome graph")

mapping_data_stats <- mapping_data %>%
  mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
  mutate(MapQ0 = Count * (MapQ >= 0)) %>% 
  mutate(MapQ30 = Count * (MapQ >= 30)) %>% 
  group_by(Reads, Method, Reference) %>%
  summarise(count = sum(Count), MapQ0 = sum(MapQ0), MapQ30 = sum(MapQ30)) %>%
  mutate(MapQ0_frac = MapQ0 / count, MapQ30_frac = MapQ30 / count) %>%
  gather("MapQ0_frac", "MapQ30_frac", key = "Filter", value = "Frac")

########

mapping_data_stats_main <- mapping_data_stats %>%
  filter(Reference != "Spliced personal graph/reference")

for (reads in unique(mapping_data_stats_main$Reads)) {

  mapping_data_stats_main_reads <- mapping_data_stats_main %>%
    filter(Reads == reads)

  mapping_data_stats_main_reads <- mapping_data_stats_main_reads %>%
    ungroup() %>%
    add_row(Reads = reads, Method = "STAR", Reference = "Spliced pangenome graph", count = 0, MapQ0 = 0, MapQ30 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
    add_row(Reads = reads, Method = "STAR", Reference = "Spliced pangenome graph", count = 0, MapQ0 = 0, MapQ30 = 0, Filter = "MapQ30_frac", Frac = 0)

  mapping_data_stats_main_reads$Method <- factor(mapping_data_stats_main_reads$Method, levels = c("HISAT2", "STAR", "vg map", "vg mpmap"))

  mapping_data_stats_main_reads$Reference = recode_factor(mapping_data_stats_main_reads$Reference,
                                                 "Spliced reference" = "Spliced\nreference",
                                                 "Spliced pangenome graph" = "Spliced\npangenome\ngraph")

  mapping_data_stats_main_reads$Filter <- recode_factor(mapping_data_stats_main_reads$Filter,
                                                         "MapQ0_frac" = "Mapped",
                                                         "MapQ30_frac" = "MapQ >= 30")

  mapping_data_stats_main_reads$FacetCol <- "Real reads,\nprimary alignments"
  mapping_data_stats_main_reads$FacetRow <- ""

  plotMappingStatsBenchmark(mapping_data_stats_main_reads, wes_cols, paste("plots/real_stats/real_r2_stats_bar_main_", reads, sep = ""))
}

########

mapping_data_stats_personal <- mapping_data_stats %>%
  filter(Reads == "ENCSR000AED_rep1") %>%
  filter(Method == "STAR" | ((Method == "vg mpmap" | Method == "Diploid reference (STAR)" | Method == "Diploid reference (STAR, LevioSAM)") & Reference == "Spliced personal graph/reference"))
  
for (reads in unique(mapping_data_stats_personal$Reads)) {
  
  mapping_data_stats_personal_reads <- mapping_data_stats_personal %>%
    filter(Reads == reads)

  mapping_data_stats_personal_reads <- mapping_data_stats_personal_reads %>%
    ungroup() %>%
    add_row(Reads = reads, Method = "vg mpmap", Reference = "Spliced reference", count = 0, MapQ0 = 0, MapQ30 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
    add_row(Reads = reads, Method = "vg mpmap", Reference = "Spliced reference", count = 0, MapQ0 = 0, MapQ30 = 0, Filter = "MapQ30_frac", Frac = 0)
  
  mapping_data_stats_personal_reads$Method <- factor(mapping_data_stats_personal_reads$Method, levels = c("STAR", "vg mpmap", "Diploid reference (STAR)", "Diploid reference (STAR, LevioSAM)"))
  
  mapping_data_stats_personal_reads$Reference = recode_factor(mapping_data_stats_personal_reads$Reference,
                                                          "Spliced reference" = "Spliced\nreference",
                                                          "Spliced personal graph/reference" = "Spliced personal\ngraph/reference")

  mapping_data_stats_personal_reads$Filter <- recode_factor(mapping_data_stats_personal_reads$Filter, 
                                                        "MapQ0_frac" = "Mapped", 
                                                        "MapQ30_frac" = "MapQ >= 30")
  
  mapping_data_stats_personal_reads$FacetCol <- "Real reads,\nprimary alignments"
  mapping_data_stats_personal_reads$FacetRow <- ""
  
  plotMappingStatsBenchmarkWide(mapping_data_stats_personal_reads, wes_cols[c(2,4,5,6)], paste("plots/real_stats/real_r2_stats_bar_personal_", reads, sep = ""))
}

########
