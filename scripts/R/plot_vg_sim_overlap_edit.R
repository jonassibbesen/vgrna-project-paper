
# plot_vg_sim_overlap_edit.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_final/")

########

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table(filename)
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Reference = dir_split[9]) 
  
  return(data)
}

overlap_threshold <- 90

########

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r3/")

overlap_data_raw_h1 <- map_dfr(list.files(path = "./methods", pattern=".*sim_vg.*_bam_ovl3_h1.txt.gz", full.names = T, recursive = T), parse_file) %>%
  select(-VarSubBP2, -VarIndelBP2) %>%
  rename(VarSubBP = VarSubBP1) %>%
  rename(VarIndelBP = VarIndelBP1) %>%
  group_by(TruthAlignmentLength, IsMapped, MapQ, GroupMapQ, Overlap, VarSubBP, VarIndelBP,	EditSubBP, EditIndelBP, PrimaryMapq, PrimaryOverlap, Type, Reads, Method, Reference) %>%
  summarise(Count = sum(Count))

overlap_data_raw_h2 <- map_dfr(list.files(path = "./methods", pattern=".*sim_vg.*_bam_ovl3_h2.txt.gz", full.names = T, recursive = T), parse_file)  %>%
  select(-VarSubBP1, -VarIndelBP1) %>%
  rename(VarSubBP = VarSubBP2) %>%
  rename(VarIndelBP = VarIndelBP2) %>%
  group_by(TruthAlignmentLength, IsMapped, MapQ, GroupMapQ, Overlap, VarSubBP, VarIndelBP,	EditSubBP, EditIndelBP, PrimaryMapq, PrimaryOverlap, Type, Reads, Method, Reference) %>%
  summarise(Count = sum(Count))

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_final/")

overlap_data <- rbind(overlap_data_raw_h1, overlap_data_raw_h2) %>%
  filter(Type == "polya_rna") %>%
  filter(Method != "star_wasp") %>%
  filter(Reference != "1kg_NA12878_gencode100") %>%
  filter(Reference != "1kg_NA12878_exons_gencode100") %>%
  filter(TruthAlignmentLength > 50) 

gc()

overlap_data %>% group_by(Reads, Reference, Method) %>%
  summarise(Count = sum(Count)) %>%
  print(n = 100)
  
overlap_data_prim1 <- overlap_data %>%
  filter(Method == "hisat2_multi10" | Method == "star_multi10" | Method == "star_wasp" | Method == "mpmap_multi10") %>%
  mutate(MapQ = PrimaryMapq) %>%
  mutate(Overlap = PrimaryOverlap) %>%
  mutate(Eval = "unique alignments")

overlap_data_prim2 <- overlap_data %>%
  filter(Method == "map_fast") %>%
  mutate(Eval = "unique alignments")

overlap_data_multi <- overlap_data %>%
  filter(Method == "hisat2_multi10" | Method == "star_multi10" | Method == "map_fast_multi10" | Method == "mpmap_multi10") %>%
  mutate(MapQ = GroupMapQ) %>%
  mutate(Eval = "multi alignments")

overlap_data <- rbind(overlap_data_prim1, overlap_data_prim2, overlap_data_multi) %>%
  mutate(Correct = Overlap >= (overlap_threshold / 100)) 
  
overlap_data$Method <- recode_factor(overlap_data$Method, 
                                     "hisat2" = "HISAT2",
                                     "hisat2_multi10" = "HISAT2",
                                     "star" = "STAR",
                                     "star_multi10" = "STAR",
                                     "map_fast" = "vg map",
                                     "map_fast_multi10" = "vg map",
                                     "mpmap" = "vg mpmap",
                                     "mpmap_multi10" = "vg mpmap")

overlap_data$Reference = recode_factor(overlap_data$Reference, 
                                         "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode80" = "80% spliced graph/reference",
                                         "gencode80" = "80% spliced graph/reference")

overlap_data$FacetCol <- paste("Simulated reads (vg), ", overlap_data$Eval, sep = "")
overlap_data$FacetRow <- ""
overlap_data$LineType <- overlap_data$Reference

overlap_data$FacetCol <- factor(overlap_data$FacetCol, levels = c("Simulated reads (vg), unique alignments", "Simulated reads (vg), multi alignments"))

gc()

########

overlap_data_edit_var <- overlap_data %>%
  mutate(Edit = VarSubBP + VarIndelBP) %>%
  mutate(Edit = ifelse(Edit > 4, 4, Edit))

for (reads in unique(overlap_data_edit_var$Reads)) {

  overlap_data_edit_var_reads <- overlap_data_edit_var %>%
    filter(Reads == reads)

  plotRocBenchmarkEdit(overlap_data_edit_var_reads, wes_cols, "Reference", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_var_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
  plotCurveBenchmarkEdit(overlap_data_edit_var_reads, wes_cols, "Reference", "Genomic variant edit distance", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_var_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
}

gc()

overlap_data_edit_error <- overlap_data %>%
  mutate(Edit = EditSubBP + EditIndelBP) %>%
  mutate(Edit = ifelse(Edit > 4, 4, Edit))

for (reads in unique(overlap_data_edit_error$Reads)) {
  
  overlap_data_edit_error_reads <- overlap_data_edit_error %>%
    filter(Reads == reads)
  
  plotRocBenchmarkEdit(overlap_data_edit_error_reads, wes_cols, "Reference", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_error_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
  plotCurveBenchmarkEdit(overlap_data_edit_error_reads, wes_cols, "Reference", "Sequencing error edit distance", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_error_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
}

gc()

overlap_data_edit_total <- overlap_data %>%
  mutate(Edit = VarSubBP + VarIndelBP + EditSubBP + EditIndelBP) %>%
  mutate(Edit = ifelse(Edit > 4, 4, Edit)) 

for (reads in unique(overlap_data_edit_total$Reads)) {
  
  overlap_data_edit_total_reads <- overlap_data_edit_total %>%
    filter(Reads == reads)
  
  plotRocBenchmarkEdit(overlap_data_edit_total_reads, wes_cols, "Reference", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_total_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
  plotCurveBenchmarkEdit(overlap_data_edit_total_reads, wes_cols, "Reference", "Edit distance", paste("plots/sim_overlap_edit/vg_sim_r3_overlap_edit_total_ovl", overlap_threshold, "_", reads, "_final", sep = ""))
}

########
