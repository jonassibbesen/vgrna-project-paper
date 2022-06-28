
# plot_vg_sim_overlap.R

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
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Reference = dir_split[9])
  
  if (grepl("_ovl0_", basename(filename))) {
    
    data <- data %>%
      add_column(Filter = "Unfiltered")
  
  } else if (grepl("_ovl3_", basename(filename))) {
    
    data <- data %>%
      add_column(Filter = "Low quality bases filtered")
    
  } else {
    
    stopifnot(FALSE)
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

overlap_threshold <- 90

########

overlap_data_raw_h1 <- map_dfr(list.files(path = "./methods", pattern=".*_bam_ovl3_h1.txt.gz", full.names = T, recursive = T), parse_file) %>%
  select(-SubstitutionBP2, -IndelBP2) %>%
  rename(SubstitutionBP = SubstitutionBP1) %>%
  rename(IndelBP = IndelBP1)

overlap_data_raw_h2 <- map_dfr(list.files(path = "./methods", pattern=".*_bam_ovl3_h2.txt.gz", full.names = T, recursive = T), parse_file)  %>%
  select(-SubstitutionBP1, -IndelBP1) %>%
  rename(SubstitutionBP = SubstitutionBP2) %>%
  rename(IndelBP = IndelBP2)

overlap_data <- rbind(overlap_data_raw_h1, overlap_data_raw_h2) %>%
  filter(Type == "polya_rna") %>%
  filter(Filter == "Low quality bases filtered") %>%
  filter(TruthAlignmentLength > 50) 

overlap_data %>% group_by(Reads, Reference, Method) %>%
  summarise(Count = sum(Count)) %>%
  print(n = 100)
  
overlap_data_prim1 <- overlap_data %>%
  filter(Method == "hisat2_multi10" | Method == "star_multi10" | Method == "star_wasp" | Method == "star_alleleseq" | Method == "star_alleleseq_levio" | Method == "mpmap_multi10") %>%
  mutate(MapQ = PrimaryMapq) %>%
  mutate(Overlap = PrimaryOverlap) %>%
  mutate(Eval = "Primary")

overlap_data_prim2 <- overlap_data %>%
  filter(Method == "map_fast") %>%
  mutate(Eval = "Primary")

overlap_data_multi <- overlap_data %>%
  filter(Method == "hisat2_multi10" | Method == "star_multi10" | Method == "map_fast_multi10" | Method == "mpmap_multi10") %>%
  mutate(MapQ = GroupMapQ) %>%
  mutate(Eval = "Mulit-mapping")

overlap_data <- rbind(overlap_data_prim1, overlap_data_prim2, overlap_data_multi) %>%
  mutate(Correct = Overlap >= (overlap_threshold / 100)) 
  
overlap_data$Method <- recode_factor(overlap_data$Method, 
                                     "hisat2" = "HISAT2",
                                     "hisat2_multi10" = "HISAT2",
                                     "star" = "STAR",
                                     "star_multi10" = "STAR",
                                     "map" = "vg map (def)", 
                                     "map_fast" = "vg map",
                                     "map_fast_multi10" = "vg map",
                                     "mpmap" = "vg mpmap",
                                     "mpmap_multi10" = "vg mpmap",
                                     "star_alleleseq" = "Diploid reference (STAR)",
                                     "star_alleleseq_levio" = "Diploid reference (STAR, LevioSAM)",
                                     "star_wasp" = "WASP (STAR)")

overlap_data[overlap_data$Method == "WASP (STAR)",]$Reference <- "1kg_NA12878_gencode100"

overlap_data$Reference = recode_factor(overlap_data$Reference, 
                                         "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                         "1kg_NA12878_gencode100" = "Spliced personal graph/reference",
                                         "1kg_NA12878_exons_gencode100" = "Spliced personal graph/reference",
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode80" = "80% spliced graph/reference",
                                         "gencode80" = "80% spliced graph/reference")

overlap_data$FacetCol <- overlap_data$Simulation
overlap_data$FacetRow <- ""
overlap_data$LineType <- overlap_data$Reference
  
########

overlap_data_debug <- overlap_data
overlap_data_debug$FacetCol <- overlap_data$Eval

for (reads in unique(overlap_data_debug$Reads)) {

  overlap_data_debug_reads <- overlap_data_debug %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQDebug(overlap_data_debug_reads, wes_cols, "Reference", paste("plots/sim_overlap/debug/vg_sim_r2_overlap_debug_ovl", overlap_threshold, "_", reads, sep = ""))
}

overlap_data_debug <- overlap_data_debug %>%
  filter(Reference == "80% spliced graph/reference") %>%
  mutate(NovelSJ = "No novel SJ") %>%
  mutate(NovelSJ = ifelse(NonAnnoSJ > 0, "Has novel SJ", NovelSJ))

overlap_data_debug$FacetRow <- overlap_data_debug$NovelSJ

for (reads in unique(overlap_data_debug$Reads)) {

  overlap_data_debug_reads <- overlap_data_debug %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQDebug(overlap_data_debug_reads, wes_cols, "Reference", paste("plots/sim_overlap/debug/vg_sim_r2_overlap_debug_ovl", overlap_threshold, "_", reads, "_novelSJ", sep = ""))
}

########

overlap_data_main <- overlap_data %>%
  filter(Method != "WASP (STAR)") %>%
  filter(Method != "Diploid reference (STAR)") %>%
  filter(Method != "vg map (def)")

overlap_data_main_prim <- overlap_data_main %>%
  filter(Eval == "Primary") %>%
  filter(Reference != "Spliced personal graph/reference") %>%
  mutate(FacetCol = paste(FacetCol, ", primary alignments", sep = ""))

for (reads in unique(overlap_data_main_prim$Reads)) {

  overlap_data_main_prim_reads <- overlap_data_main_prim %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_prim_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_ovl", overlap_threshold, "_", reads, "_primary", sep = ""))
}

overlap_data_main_multi <- overlap_data_main %>%
  filter(Eval == "Mulit-mapping") %>%
  filter(Reference != "Spliced personal graph/reference") %>%
  mutate(FacetCol = paste(FacetCol, ", multi alignments", sep = ""))

for (reads in unique(overlap_data_main_multi$Reads)) {

  overlap_data_main_multi_reads <- overlap_data_main_multi %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

########

overlap_data_main_multi_novel_sj <- overlap_data_main_multi %>%
  filter(Reference == "80% spliced graph/reference") %>%
  mutate(NovelSJ = "No novel SJ") %>%
  mutate(NovelSJ = ifelse(NonAnnoSJ > 0, "Has novel SJ", NovelSJ)) %>%
  mutate(LineType = NovelSJ)

overlap_data_main_multi_novel_sj$LineType <- recode_factor(overlap_data_main_multi_novel_sj$LineType,
                                     "No novel SJ" = "No novel splice-junctions     ",
                                     "Has novel SJ" = "Has novel splice-junctions    ")

for (reads in unique(overlap_data_main_multi_novel_sj$Reads)) {

  overlap_data_main_multi_novel_sj_reads <- overlap_data_main_multi_novel_sj %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_novel_sj_reads, wes_cols, "Reads", paste("plots/sim_overlap/vg_sim_r2_overlap_main_ovl", overlap_threshold, "_", reads, "_multi_novel_sj", sep = ""))
}

########

overlap_data_main <- overlap_data_main %>%
  filter(Reads == "sim_vg_r2_ENCSR000AED_rep1_uni")

overlap_data_main_error <- overlap_data_main %>%
  filter(Eval == "Primary") %>%
  filter(Reference != "Spliced personal graph/reference") %>%
  filter(Reference != "80% spliced graph/reference") %>%
  filter(Reference != "Spliced reference" | Method == "STAR") %>%
  mutate(FacetCol = paste(FacetCol, ", primary alignments", sep = ""))

overlap_data_main_error90 <- overlap_data_main_error %>%
  mutate(FacetCol = paste(FacetCol, "\nTrue alignment cover >= 90%", sep = ""))

for (reads in unique(overlap_data_main_error90$Reads)) {

  overlap_data_main_error90_reads <- overlap_data_main_error90 %>%
    filter(Reads == reads)

  plotErrorBenchmark(overlap_data_main_error90_reads, wes_cols, paste("plots/sim_overlap/vg_sim_r2_overlap_main_ovl", overlap_threshold, "_", reads, "_primary", sep = ""))
}

overlap_data_main_error1 <- overlap_data_main_error %>%
  mutate(FacetCol = paste(FacetCol, "\nTrue alignment cover > 0%", sep = "")) %>%
  mutate(Correct = Overlap > 0)

for (reads in unique(overlap_data_main_error1$Reads)) {

  overlap_data_main_error1_reads <- overlap_data_main_error1 %>%
    filter(Reads == reads)

  plotErrorBenchmark(overlap_data_main_error1_reads, wes_cols, paste("plots/sim_overlap/vg_sim_r2_overlap_main_ovl1_", reads, "_primary", sep = ""))
}

########

addReadsFracToFacetCol <- function(data, total_num_reads) {

  num_reads <- as.numeric(unique(data %>%
    group_by(Reads, Method, Reference) %>%
    summarise(Count = sum(Count)) %>%
    ungroup() %>%
    select(Count)))

  data <- data %>%
    mutate(FacetCol = paste(FacetCol, ", ", signif(num_reads / total_num_reads * 100, 3), "% reads", sep = ""))

  return(data)
}

overlap_data_main_multi <- overlap_data_main %>%
  filter(Eval == "Mulit-mapping") %>%
  filter(Reference != "Spliced personal graph/reference") %>%
  mutate(FacetCol = paste(FacetCol, ", multi alignments", sep = ""))

total_num_reads <- as.numeric(unique(overlap_data_main_multi %>%
  group_by(Reads, Method, Reference) %>%
  summarise(Count = sum(Count)) %>%
  ungroup() %>%
  select(Count)))

overlap_data_main_multi_nov <- overlap_data_main_multi %>%
  filter(SubstitutionBP == 0 & IndelBP == 0) %>%
  mutate(FacetCol = paste(FacetCol, "\nNo variants", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_nov$Reads)) {

  overlap_data_main_multi_nov_reads <- overlap_data_main_multi_nov %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_nov_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_nov_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

overlap_data_main_multi_snv1 <- overlap_data_main_multi %>%
  filter(SubstitutionBP == 1 & IndelBP == 0) %>%
  mutate(FacetCol = paste(FacetCol, "\n1 SNV (no indels)", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_snv1$Reads)) {

  overlap_data_main_multi_snv1_reads <- overlap_data_main_multi_snv1 %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_snv1_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_snv1_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

overlap_data_main_multi_snv2 <- overlap_data_main_multi %>%
  filter(SubstitutionBP == 2 & IndelBP == 0) %>%
  mutate(FacetCol = paste(FacetCol, "\n2 SNVs (no indels)", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_snv2$Reads)) {

  overlap_data_main_multi_snv2_reads <- overlap_data_main_multi_snv2 %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_snv2_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_snv2_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

overlap_data_main_multi_snv3 <- overlap_data_main_multi %>%
  filter(SubstitutionBP == 3 & IndelBP == 0) %>%
  mutate(FacetCol = paste(FacetCol, "\n3 SNVs (no indels)", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_snv3$Reads)) {

  overlap_data_main_multi_snv3_reads <- overlap_data_main_multi_snv3 %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_snv3_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_snv3_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

overlap_data_main_multi_snv4 <- overlap_data_main_multi %>%
  filter(SubstitutionBP > 3 & IndelBP == 0) %>%
  mutate(FacetCol = paste(FacetCol, "\n>3 SNVs (no indels)", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_snv4$Reads)) {

  overlap_data_main_multi_snv4_reads <- overlap_data_main_multi_snv4 %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_snv4_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_snv4_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

overlap_data_main_multi_indel <- overlap_data_main_multi %>%
  filter(IndelBP > 0) %>%
  mutate(FacetCol = paste(FacetCol, "\n>0 indels", sep = "")) %>%
  addReadsFracToFacetCol(total_num_reads)

for (reads in unique(overlap_data_main_multi_indel$Reads)) {

  overlap_data_main_multi_indel_reads <- overlap_data_main_multi_indel %>%
    filter(Reads == reads)

  plotRocBenchmarkMapQ(overlap_data_main_multi_indel_reads, wes_cols, "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_main_indel_ovl", overlap_threshold, "_", reads, "_multi", sep = ""))
}

########

overlap_data_personal_prim <- overlap_data %>%
  filter(Reads == "sim_vg_r2_ENCSR000AED_rep1_uni") %>%
  filter((Method == "STAR" & Reference == "Spliced reference") | ((Method == "vg mpmap" | Method == "Diploid reference (STAR)" | Method == "Diploid reference (STAR, LevioSAM)") & Reference == "Spliced personal graph/reference")) %>%
  filter(Eval == "Primary") %>%
  mutate(FacetCol = paste(FacetCol, ", primary alignments", sep = ""))

for (reads in unique(overlap_data_personal_prim$Reads)) {
  
  overlap_data_personal_prim_reads <- overlap_data_personal_prim %>%
    filter(Reads == reads)
  
  plotRocBenchmarkMapQ(overlap_data_personal_prim_reads, wes_cols[c(2,4,5,6)], "Reference", paste("plots/sim_overlap/vg_sim_r2_overlap_personal_ovl", overlap_threshold, "_", reads, "_primary", sep = ""))
}

########
