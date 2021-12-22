
# parse_expression_data.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")
library("truncnorm")

# source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########

dataset <- "SRR1153470"
sim_mean <- 277
sim_sd <- 43

# dataset <- "ENCSR000AED_rep1"
# sim_mean <- 216
# sim_sd <- 24

#read_type <- "sim_vg"
read_type <- "real"

#ref_name <- "1kg_NA12878_gencode100"
#ref_name <- "1kg_EURnonCEU_af002_gencode100"
#ref_name <- "1kg_EURnonCEU_af002_gencode100_unidi"
ref_name <- "1kg_nonCEU_af001_gencode100"
#ref_name <- "1kg_nonCEU_af001_gencode100_unidi"
#ref_name <- "1kg_all_af001_gencode100"
#ref_name <- "1kg_all_af001_gencode100_unidi"

hap_prob_thres <- 0.8

update_tpm <- function(exp_data) {
  
  exp_data <- exp_data %>%
    mutate(TPM = count / (length - etruncnorm(1, length, sim_mean, sim_sd))) %>%
    replace_na(list(TPM = 0)) %>%
    mutate(TPM = 10^6 * TPM / sum(TPM))
  
  return(exp_data)
}

parse_salmon <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = Name, tpm_est = TPM, count_est = NumReads, length = Length) %>%
    select(-EffectiveLength)
  
  if ("FracNonZero" %in% names(data)) {
    
    data <- data %>%
      select(-FracNonZero)     
  }
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-eff_length)
  
  return(data)
}

parse_rsem <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = transcript_id, tpm_est = TPM, count_est = expected_count) %>%
    select(-gene_id, -effective_length, -FPKM, -IsoPct)
  
  if ("posterior_mean_count" %in% names(data)) {
    
    data <- data %>%
      mutate(count_est = posterior_mean_count) %>%
      mutate(tpm_est = pme_TPM) %>%
      select(-posterior_mean_count, -posterior_standard_deviation_of_count, -pme_TPM, -pme_FPKM, -IsoPct_from_pme_TPM)     
  }
  
  return(data)
}

parse_rpvg <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount, length = Length) %>%
    select(-EffectiveLength)
  
  if ("ClusterRelativeExpression" %in% names(data)) {
    
    data <- data %>%
      select(-ClusterRelativeExpression)     
  }
  
  if (grepl("_multi_", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "multi", sep = "_"))
  } 
  
  return(data)
}

getStats <- function(data) {

  data_stats <- data %>%
    mutate(ard_count = abs(count_est - count_sim) / (count_est + count_sim)) %>%
    mutate(ard_tpm = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
    replace_na(list(ard_count = 0)) %>%
    replace_na(list(ard_tpm = 0)) %>%
    group_by(Reads, Method, Graph) %>%
    summarise(
      n = n(), 
      Expressed = sum(tpm_est > 0), 
      ExpCorrect = mean((tpm_est > 0) == (tpm_sim > 0)), 
      ExpSimCorrect = sum((tpm_est > 0) & (tpm_sim > 0)) / sum(tpm_sim > 0), 
      Pearson_count = cor(count_est, count_sim, method = "pearson"),
      Pearson_tpm = cor(tpm_est, tpm_sim, method = "pearson"),
      Spearman_count = cor(count_est, count_sim, method = "spearman"),
      Spearman_tpm = cor(tpm_est, tpm_sim, method = "spearman"),
      ARD_mean_count = mean(ard_count),
      ARD_mean_tpm = mean(ard_tpm),
      frac_hap_error_count = sum((!is_hap) * count_est) / sum(count_est),
      frac_hap_error_tpm = sum((!is_hap) * tpm_est) / sum(tpm_est),
      num_hap_error_count = sum((!is_hap) & count_est > 0),
      num_hap_error_tpm = sum((!is_hap) & tpm_est > 0))
  
  return(data_stats)
}

# sim_exp_h1 <- read_table2(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h1_paths.txt.gz", sep = "")) %>%
#   group_by(path) %>%
#   summarise(count = n()) %>%
#   mutate(count = count / 2)
# 
# sim_exp_h2 <- read_table2(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h2_paths.txt.gz", sep = "")) %>%
#   group_by(path) %>%
#   summarise(count = n()) %>%
#   mutate(count = count / 2)
# 
# sim_exp <- bind_rows(sim_exp_h1, sim_exp_h2)
# save(sim_exp, file = paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))


identical_seqs <- read_table2(paste("../quant/graphs/1kg_NA12878_exons_gencode100_allpaths/", ref_name, "_hst_overlap.txt", sep = ""))
rsem <- read_table2(paste("../quant/rsem/", dataset, "/1kg_NA12878_gencode100_", dataset , "_rsem.isoforms.results", sep = ""))

load(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))

sim_exp <- sim_exp %>%
  right_join(rsem, by = c("path" = "transcript_id")) %>%
  left_join(identical_seqs, by = c("path" = "Name1")) %>%
  rename(name = Name2) %>%
  mutate(name = ifelse(is.na(name), paste(path, "_na", sep = ""), name)) %>%
  replace_na(list(count = 0)) %>%
  update_tpm() %>%
  group_by(name) %>%
  summarise(length_sim = max(length), tpm_sim = sum(TPM), count_sim = sum(count))

sim_exp <- sim_exp %>%
  ungroup() %>%
  mutate(count_sim = count_sim / sum(count_sim) * 1000000)
  
if (read_type == "real") {
  
  sim_exp <- sim_exp %>%
    mutate(tpm_sim = 1) %>%
    mutate(count_sim = 1)
}

files <- c(list.files(path = "methods", pattern = "rpvg.*.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "quant.sf.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "abundance.tsv.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "isoforms.results.gz", full.names = T, recursive = T))

for (f in files) { 

  if (grepl("_joint.txt.gz", f)) {
    
    next 
  }
  
  if (grepl("_gibbs.txt.gz", f)) {
    
    next 
  }
  
  if (!grepl(dataset, f)) {
    
    next 
  }
  
  if (!grepl(read_type, f)) {
    
    next 
  }
   
  if (!grepl(ref_name, f)) {
    
    next 
  }
  
  if (grepl("unidi", f) != grepl("unidi", ref_name)) {
    
    next 
  }
  
  print(f)
  
  if (grepl("quant.sf", f) | grepl("quant_boot.sf", f)) {
    
    exp_data <- parse_salmon(f)
  
  } else if (grepl("abundance.tsv", f)) {

    exp_data <- parse_kallisto(f)

  } else if (grepl("isoforms.results", f)) {

    exp_data <- parse_rsem(f)
    
  } else {
    
    exp_data <- parse_rpvg(f)
  }

  exp_data <- exp_data %>% 
    ungroup() %>%
    mutate(count_est = count_est / sum(count_est) * 1000000)
  
  # exp_data %>% filter(tpm_est > 0) %>% group_by(Reads, Method, Graph) %>% summarise(sum_count = sum(count_est), min_count = min(count_est), sum_tpm = sum(tpm_est), min_tpm = min(tpm_est)) %>% print(n = 100)
  
  if (read_type == "sim_vg") {
    
    exp_data <- exp_data %>%
      rename(count = count_est, TPM = tpm_est) %>%
      update_tpm() %>%
      rename(count_est = count, tpm_est = TPM)
  }
  
  # sim_exp %>% filter(tpm_sim > 0) %>% summarise(sum_count = sum(count_sim), min_count = min(count_sim), sum_tpm = sum(tpm_sim), min_tpm = min(tpm_sim)) %>% print(n = 100)
  # exp_data %>% filter(tpm_est > 0) %>% group_by(Reads, Method, Graph) %>% summarise(sum_count = sum(count_est), min_count = min(count_est), sum_tpm = sum(tpm_est), min_tpm = min(tpm_est)) %>% print(n = 100)
  
  exp_data <- exp_data %>% 
    full_join(sim_exp, by = "name") %>%
    mutate(is_hap = ifelse(is.na(tpm_sim), F, T)) %>%
    replace_na(list(HaplotypeProbability = 0, tpm_est = 0, count_est = 0, Reads = exp_data$Reads[1], Method = exp_data$Method[1], Graph = exp_data$Graph[1], tpm_sim = 0, count_sim = 0)) %>%
    separate(name, c("transcript", "hap_id"), "_") 
  
  
  # exp_data_debug <- exp_data %>% group_by(ClusterID) %>% mutate(ClusterSize = n()) %>% filter(HaplotypeProbability >= 0.8 | is_hap) 
  # exp_data_debug <- exp_data_debug %>% group_by(transcript) %>% mutate(num_hap = sum(is_hap), max_hap_len = max(is_hap * length))
  # 
  # exp_data_debug %>% filter(tpm_est > 0) %>% filter(!is_hap) %>% filter(num_hap == 1) %>% filter(abs(length - max_hap_len) == 1) %>% arrange(desc(tpm_est)) %>% print(n = 70)
  # 
  # exp_data_debug %>% filter(transcript == "ENST00000646664.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000514057.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000253788.11") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000380394.8") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000436459.2") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000414273.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000457540.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000417615.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000471152.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000309311.6") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000525807.5") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000393820.2") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000322723.8") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000330899.4") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000360004.5") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000300026.3") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000374975.3") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000557016.5") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000330459.7") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000234313.7") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000366560.3") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # 
  # exp_data_debug %>% filter(tpm_est > 0) %>% filter(!is_hap) %>% filter(num_hap == 1) %>% filter(abs(length - max_hap_len) == 1) %>% arrange(desc(tpm_est)) %>% print(n = 20)
  # 
  # exp_data_debug %>% filter(transcript == "ENST00000335895.12") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000395839.5") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000395837.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000505490.2") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000252725.10") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000322428.9") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000264156.2") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000297258.10") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000432629.1") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  # exp_data_debug %>% filter(transcript == "ENST00000260563.4") %>% filter(is_hap | (!is_hap & tpm_est > 0))
  
  
  exp_data_hap_prob <- exp_data %>%
    group_by(HaplotypeProbability, Reads, Method, Graph) %>%
    summarise(TP = sum((tpm_sim > 0) & (tpm_est > 0)),
              TN = sum((tpm_sim == 0) & (tpm_est == 0)),
              FP = sum((tpm_sim == 0) & (tpm_est > 0)),
              FN = sum((tpm_sim > 0) & (tpm_est == 0)),
              TP_tpm = sum((tpm_sim > 0) * tpm_est),
              FP_tpm = sum((tpm_sim == 0) * tpm_est))
  
  exp_data <- exp_data %>% 
    mutate(count_est = ifelse(HaplotypeProbability >= hap_prob_thres, count_est, 0)) %>%
    mutate(tpm_est = ifelse(HaplotypeProbability >= hap_prob_thres, tpm_est, 0))
  
  exp_data_hap_exp <- exp_data %>%
    group_by(tpm_est, Reads, Method, Graph) %>%
    summarise(TP = sum((tpm_sim > 0) & (tpm_est > 0)),
              TN = sum((tpm_sim == 0) & (tpm_est == 0)),
              FP = sum((tpm_sim == 0) & (tpm_est > 0)),
              FN = sum((tpm_sim > 0) & (tpm_est == 0)),
              TP_tpm = sum((tpm_sim > 0) * tpm_est),
              FP_tpm = sum((tpm_sim == 0) * tpm_est))
  
  exp_data_stats_all <- exp_data %>%
    getStats() %>%
    add_column(Type = "All")

  exp_data_stats_hap <- exp_data %>%
    filter(is_hap) %>%
    getStats() %>%
    add_column(Type = "Haplotype")
  
  exp_data_stats_txp <- exp_data %>%
    group_by(transcript, Reads, Method, Graph) %>%
    summarise(count_est = sum(count_est), count_sim = sum(count_sim), tpm_est = sum(tpm_est), tpm_sim = sum(tpm_sim), is_hap = max(is_hap)) %>%
    ungroup() %>%
    getStats() %>%
    add_column(Type = "Transcript")
  
  exp_data_hap_exp <- exp_data_hap_exp %>%
    add_column(hap_prob_thres = hap_prob_thres) 
  
  exp_data_stats <- rbind(exp_data_stats_all, exp_data_stats_hap, exp_data_stats_txp) %>%
    add_column(hap_prob_thres = hap_prob_thres) 
  
  print(exp_data_stats)
  
  save(exp_data_hap_prob, exp_data_hap_exp, exp_data_stats, file = paste("rdata/", exp_data_stats$Method[1], exp_data_stats$Reads[1], exp_data_stats$Graph[1], ".RData", sep = "", collapse = ""))
}
