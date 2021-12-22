
# plot_expression_correlation_rsem.R

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


parse_rpvg <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount) %>%
    # mutate(count_est = ifelse(HaplotypePosterior >= 0, count_est, 0)) %>%
    # mutate(tpm_est = count_est / as.double(EffectiveLength)) %>%
    # replace_na(list(tpm_est = 0)) %>%
    # mutate(tpm_est = 10^6 * tpm_est / sum(tpm_est)) %>%
    mutate(count_est = ifelse(ClusterRelativeExpression < 10^-8, 0, count_est)) %>%
    mutate(tpm_est = ifelse(ClusterRelativeExpression < 10^-8, 0, tpm_est)) %>%
    select(-ClusterID, -Length, -EffectiveLength, -ClusterRelativeExpression)
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    add_column(HaplotypePosterior = 1) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-length, -eff_length)
  
  return(data)
}

parse_salmon <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    add_column(HaplotypePosterior = 1) %>%
    rename(name = Name, tpm_est = TPM, count_est = NumReads) %>%
    select(-Length, -EffectiveLength)
  
  return(data)
}

identical_seqs <- read_table2("graphs/1kg_NA12878_exons_gencode100_allpaths/1kg_all_af001_gencode100_hst_overlap.txt")

sim_exp_h1 <- read_table2("sim_1kg_NA12878_gencode100/SRR1153470/rsem/sim_1kg_NA12878_gencode100_SRR1153470_rsem_h1.sim.isoforms.results")
sim_exp_h2 <- read_table2("sim_1kg_NA12878_gencode100/SRR1153470/rsem/sim_1kg_NA12878_gencode100_SRR1153470_rsem_h2.sim.isoforms.results")

sim_exp <- bind_rows(sim_exp_h1, sim_exp_h2) %>%
  left_join(identical_seqs, by = c("transcript_id" = "Name1")) %>%
  rename(name = Name2) %>%
  mutate(name = ifelse(is.na(name), paste(transcript_id, "_na", sep = ""), name)) %>%
  mutate(TPM = count / effective_length) %>%
  replace_na(list(TPM = 0)) %>%
  mutate(TPM = 10^6 * TPM / sum(TPM)) %>%
  group_by(name) %>%
  summarise(tpm_sim = sum(TPM), count_sim = sum(count))

#for (f in list.files(pattern = ".*_map_1kg_all_af001_gencode100.*sim_rsem_SRR1153470.txt.gz", full.names = T, recursive = T)) { 

#for (f in list.files(path = "./methods/kallisto/expression/polya_rna/sim_rsem/SRR1153470/kallisto/1kg_all_af001_gencode100/", pattern = ".*abundance.tsv.gz", full.names = T, recursive = T)) { 
  
for (f in list.files(path = "./methods/salmon/expression/polya_rna/sim_rsem/SRR1153470/salmon/1kg_all_af001_gencode100/", pattern = ".*quant.sf.gz", full.names = T, recursive = T)) { 
  
  print(f)
  exp_data <- parse_salmon(f)
  
  exp_data <- exp_data %>%
    add_column(Type = "sim_rsem")
  
  #exp_data <- map_dfr(list.files(pattern = ".*rpvg_exact_mpmap.*1kg_all_af001_gencode100_genes_sim.*.txt", full.names = T, recursive = T), parse_rpvg)
  #exp_data <- map_dfr(list.files(path = "./methods/kallisto/expression/polya_rna/sim_rsem/SRR1153470/kallisto/1kg_nonCEU_af001_gencode100/", pattern = ".*abundance.tsv", full.names = T, recursive = T), parse_kallisto)
  #exp_data <- map_dfr(list.files(path = "./methods/salmon/expression/polya_rna/sim_rsem/SRR1153470/salmon/1kg_NA12878_gencode100/", pattern = ".*quant.sf", full.names = T, recursive = T), parse_salmon)
  
  exp_data_sim <- exp_data %>% 
    full_join(sim_exp, by = "name") %>%
    mutate(is_hap = ifelse(is.na(tpm_sim), F, T)) %>%
    replace_na(list(HaplotypePosterior = 0, count_est = 0, tpm_est = 0, Type = exp_data$Type[1], Reads = exp_data$Reads[1], Method = exp_data$Method[1], Graph = exp_data$Graph[1], tpm_sim = 0, count_sim = 0))
  
  print(nrow(exp_data_sim))
  
  exp_data_sim %>%
    filter(is_hap == FALSE) %>%
    filter(tpm_est > 0)  %>%
    select(-count_sim, -tpm_sim, -Graph, -Reads, -Method, -Type) %>%
    mutate(tpm_est_frac = tpm_est / sum(tpm_est)) %>%
    arrange(desc(tpm_est)) %>%
    mutate(tpm_est_cs_frac = cumsum(tpm_est_frac)) %>%
    print(n = 10)
  
  exp_data_hap_pos <- exp_data_sim %>%
    group_by(HaplotypePosterior, Type, Reads, Method, Graph) %>%
    summarise(TP_hap = sum(is_hap), TP_hap_tpm = sum(tpm_est * is_hap), 
              FP_hap = sum(!is_hap), FP_hap_tpm = sum(tpm_est * !is_hap), 
              TP_sim = sum(tpm_sim > 0), TP_sim_tpm = sum(tpm_est * (tpm_sim > 0)), 
              FP_sim = sum(tpm_sim == 0), FP_sim_tpm = sum(tpm_est * (tpm_sim == 0)), 
              TP_pos = sum((tpm_sim > 0) & (tpm_est > 0)),
              TN_pos = sum((tpm_sim == 0) & (tpm_est == 0)),
              FP_pos = sum((tpm_sim == 0) & (tpm_est > 0)),
              FN_pos = sum((tpm_sim > 0) & (tpm_est == 0)))
  
  exp_data_sim <- exp_data_sim %>% 
    mutate(count_est = ifelse(HaplotypePosterior >= 0.9, count_est, 0)) %>%
    mutate(tpm_est = ifelse(HaplotypePosterior >= 0.9, tpm_est, 0)) %>%
    select(-name, -HaplotypePosterior, -count_est, -count_sim) 
  
  exp_data_ic_count <- exp_data_sim %>%
    filter(!is_hap & tpm_est > 0) %>%
    select(-tpm_sim, is_hap) 
  
  print(nrow(exp_data_ic_count))
  
  exp_data_all_stats <- exp_data_sim %>%
    mutate(ard = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
    replace_na(list(ard = 0)) %>%
    mutate(non_hap_tpm = tpm_est * !is_hap) %>%
    group_by(Type, Reads, Method, Graph) %>%
    summarise(n = n(), Pearson = cor(tpm_est, tpm_sim, method = "pearson"), PearsonLog = cor(log(tpm_est + 1), log(tpm_sim + 1), method = "pearson"), Spearman = cor(tpm_est, tpm_sim, method = "spearman"), ARD_mean = mean(ard), ARD_median = median(ard), hap_error = sum(non_hap_tpm) / sum(tpm_est), num_hap_error = sum(!is_hap & tpm_est > 0)) %>%
    add_column(is_hap = NA)
  
  exp_data_all_stats_hap <- exp_data_sim %>%
    mutate(ard = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
    replace_na(list(ard = 0)) %>%
    mutate(non_hap_tpm = tpm_est * !is_hap) %>%
    group_by(Type, Reads, Method, Graph, is_hap) %>%
    summarise(n = n(), Pearson = cor(tpm_est, tpm_sim, method = "pearson"), PearsonLog = cor(log(tpm_est + 1), log(tpm_sim + 1), method = "pearson"), Spearman = cor(tpm_est, tpm_sim, method = "spearman"), ARD_mean = mean(ard), ARD_median = median(ard), hap_error = sum(non_hap_tpm) / sum(tpm_est), num_hap_error = sum(!is_hap & tpm_est > 0))
  
  exp_data_all_stats <- rbind(exp_data_all_stats, exp_data_all_stats_hap)
  
  save(exp_data_hap_pos, exp_data_ic_count, exp_data_all_stats, file = paste("rdata/", exp_data_sim$Type[1], exp_data_sim$Method[1], exp_data_sim$Reads[1], exp_data_sim$Graph[1], ".RData", sep = "", collapse = ""))

}