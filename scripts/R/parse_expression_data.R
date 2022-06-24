
# parse_expression_data.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_r2/")

########

update_tpm <- function(exp_data, sim_mean, sim_sd) {
  
  exp_data <- exp_data %>%
    mutate(TPM = count / (length - etruncnorm(1, length, sim_mean, sim_sd))) %>%
    replace_na(list(TPM = 0)) %>%
    mutate(TPM = 10^6 * TPM / sum(TPM))
  
  return(exp_data)
}

parse_salmon <- function(filename) {
  
  data <- read_table(gzfile(filename), col_names = T)
  data <- data %>%
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
  
  data <- read_table(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-eff_length)
  
  return(data)
}

parse_rsem <- function(filename) {
  
  data <- read_table(gzfile(filename), col_names = T)
  data <- data %>%
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
  
  data <- read_table(gzfile(filename), col_names = T)
  data <- data %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount, length = Length) %>%
    select(-ClusterID, -EffectiveLength)
  
  if (!"HaplotypeProbability" %in% names(data)) {
    
    data <- data %>%
      add_column(HaplotypeProbability = 1)     
  }
  
  return(data)
}


getStats <- function(data) {

  data_stats <- data %>%
    mutate(ard_count = abs(count_est - count_truth) / (count_est + count_truth)) %>%
    mutate(ard_tpm = abs(tpm_est - tpm_truth) / (tpm_est + tpm_truth)) %>%
    replace_na(list(ard_count = 0)) %>%
    replace_na(list(ard_tpm = 0)) %>%
    summarise(
      n = n(), 
      Expressed = sum(tpm_est > 0), 
      ExpCorrect = mean((tpm_est > 0) == (tpm_truth > 0)), 
      ExpSimCorrect = sum((tpm_est > 0) & (tpm_truth > 0)) / sum(tpm_truth > 0), 
      Pearson_count = cor(count_est, count_truth, method = "pearson"),
      Pearson_tpm = cor(tpm_est, tpm_truth, method = "pearson"),
      Spearman_count = cor(count_est, count_truth, method = "spearman"),
      Spearman_tpm = cor(tpm_est, tpm_truth, method = "spearman"),
      ARD_mean_count = mean(ard_count),
      ARD_mean_tpm = mean(ard_tpm),
      frac_hap_error_count = sum((!is_hap) * count_est) / sum(count_est),
      frac_hap_error_tpm = sum((!is_hap) * tpm_est) / sum(tpm_est),
      num_hap_error_count = sum((!is_hap) & count_est > 0),
      num_hap_error_tpm = sum((!is_hap) & tpm_est > 0))
  
  return(data_stats)
}

hap_prob_thres <- 0.8

decoy_transcripts <- rbind(read_table("decoys/MT/gencode100_MT.txt.gz"), read_table("decoys/SCA/gencode100_SCA.txt.gz")) %>% 
  select(Transcript) %>%
  rename(name = Transcript)

parse_data <- function(dataset, sim_mean, sim_sd, read_type, ref_name) {
 
  # truth_exp_h1 <- read_table(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h1_paths.txt.gz", sep = "")) %>%
  #   group_by(path) %>%
  #   summarise(count = n()) %>%
  #   mutate(count = count / 2)
  # 
  # truth_exp_h2 <- read_table(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h2_paths.txt.gz", sep = "")) %>%
  #   group_by(path) %>%
  #   summarise(count = n()) %>%
  #   mutate(count = count / 2)
  # 
  # truth_exp <- bind_rows(truth_exp_h1, truth_exp_h2)
  # save(truth_exp, file = paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))
  
  identical_seqs <- read_table(paste("../quant_r1/graphs/1kg_NA12878_exons_gencode100_allpaths/", ref_name, "_hst_overlap.txt", sep = ""))
  rsem <- read_table(paste("../quant_r1/rsem/", dataset, "/1kg_NA12878_gencode100_", dataset , "_rsem.isoforms.results", sep = "")) %>%
    select(-effective_length, -expected_count, -TPM, -FPKM, -IsoPct)
  
  load(paste("../quant_r1/sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))
  
  truth_exp <- sim_exp %>%
    right_join(rsem, by = c("path" = "transcript_id")) %>%
    left_join(identical_seqs, by = c("path" = "Name1")) %>%
    replace_na(list(count = 0)) %>%
    group_by(path) %>% 
    mutate(n = n()) %>%
    mutate(count = count / n) %>%
    select(-n) %>%
    ungroup() %>%
    rename(name = Name2) %>%
    mutate(name = ifelse(is.na(name), paste(path, "_na", sep = ""), name)) %>%
    update_tpm(sim_mean, sim_sd) %>%
    group_by(name) %>%
    summarise(tpm_truth = sum(TPM), count_truth = sum(count))
  
  truth_exp <- truth_exp %>%
    ungroup() %>%
    mutate(count_truth = count_truth / sum(count_truth) * 1000000)
  
  if (read_type == "real") {
    
    truth_exp <- truth_exp %>%
      mutate(tpm_truth = rep(c(1,2), length.out = nrow(truth_exp))) %>%
      mutate(count_truth = rep(c(1,2), length.out = nrow(truth_exp)))
  }
  
  print(sum(truth_exp$count_truth))
  print(sum(truth_exp$tpm_truth))
  
  gc()
  
  write.table(truth_exp, file = paste("truth/truth_exp_", read_type, "_", dataset, "_", ref_name, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
  
  files <- c(list.files(path = "methods", pattern = "rpvg.*.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "quant.sf.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "abundance.tsv.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "isoforms.results.gz", full.names = T, recursive = T))
  
  for (f in files) { 
    
    if (grepl("joint", f)) {
      
      next 
    }
    
    if (grepl("gibbs", f)) {
      
      next 
    }
    
    if (!grepl(dataset, f)) {
      
      next 
    }
    
    if (!grepl(read_type, f)) {
      
      next 
    }
    
    if (!grepl(paste("/", dataset, "/", sep = ""), f)) {
      
      next 
    }
    
    if (!grepl(paste("/", read_type, sep = ""), f)) {
      
      next 
    }
    
    if (!grepl(paste("/", ref_name, sep = ""), f)) {
      
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
   
    print(sum(exp_data$count_est))
    print(sum(exp_data$tpm_est))
    
    exp_data <- exp_data %>% 
      filter(name != "Unknown")
    
    if (read_type == "sim_vg") {
      
      exp_data <- exp_data %>%
        rename(count = count_est, TPM = tpm_est) %>%
        update_tpm(sim_mean, sim_sd) %>%
        rename(count_est = count, tpm_est = TPM)
    }
    
    print(nrow(exp_data))
    print(nrow(decoy_transcripts))
    
    exp_data <- exp_data %>% 
      filter(!(name %in% decoy_transcripts$name))
    
    print(nrow(exp_data))
    print(nrow(truth_exp))
    
    exp_data <- exp_data %>% 
      select(-length) %>%
      ungroup() %>%
      mutate(count_est = count_est / sum(count_est) * 1000000) %>%
      full_join(truth_exp, by = "name") %>%
      mutate(is_hap = ifelse(is.na(tpm_truth), F, T)) %>%
      replace_na(list(HaplotypeProbability = 0, tpm_est = 0, count_est = 0, tpm_truth = 0, count_truth = 0)) %>%
      mutate(transcript = gsub("_[0-9na_]+$", "", name)) %>%
      select(-name)

    exp_data %>% 
      filter(grepl("_", transcript)) %>%
      filter(!grepl("_PAR_Y", transcript)) %>%
      print()
    
    print(nrow(exp_data))
    
    gc()
    
    exp_data_hap_prob <- exp_data %>%
      mutate(HaplotypeProbability = signif(HaplotypeProbability, 3)) %>%
      group_by(HaplotypeProbability) %>%
      summarise(TP = sum((tpm_truth > 0) & (tpm_est > 0)),
                TN = sum((tpm_truth == 0) & (tpm_est == 0)),
                FP = sum((tpm_truth == 0) & (tpm_est > 0)),
                FN = sum((tpm_truth > 0) & (tpm_est == 0)),
                TP_tpm = sum((tpm_truth > 0) * tpm_est),
                FP_tpm = sum((tpm_truth == 0) * tpm_est))
    
    print(nrow(exp_data_hap_prob))
    
    exp_data_stats_txp <- exp_data %>%
      group_by(transcript) %>%
      summarise(count_est = sum(count_est), count_truth = sum(count_truth), tpm_est = sum(tpm_est), tpm_truth = sum(tpm_truth), is_hap = max(is_hap)) %>%
      ungroup() %>%
      getStats() %>%
      add_column(Type = "Transcript")
    
    exp_data <- exp_data %>% 
      mutate(count_est = ifelse(HaplotypeProbability >= hap_prob_thres, count_est, 0)) %>%
      mutate(tpm_est = ifelse(HaplotypeProbability >= hap_prob_thres, tpm_est, 0))
    
    exp_data_stats_all <- exp_data %>%
      getStats() %>%
      add_column(Type = "All")
    
    exp_data_stats_exp <- exp_data %>%
      filter(tpm_est > 0) %>%
      getStats() %>%
      add_column(Type = "Expressed")
    
    exp_data_stats_hap <- exp_data %>%
      filter(is_hap) %>%
      getStats() %>%
      add_column(Type = "Haplotype")
    
    exp_data_stats <- rbind(exp_data_stats_all, exp_data_stats_exp, exp_data_stats_hap, exp_data_stats_txp) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
  
    print(nrow(exp_data_stats))
    
    exp_data_hap_exp <- exp_data %>%
      mutate(tpm_est = signif(tpm_est, 3)) %>%
      group_by(tpm_est) %>%
      summarise(TP = sum((tpm_truth > 0) & (tpm_est > 0)),
                TN = sum((tpm_truth == 0) & (tpm_est == 0)),
                FP = sum((tpm_truth == 0) & (tpm_est > 0)),
                FN = sum((tpm_truth > 0) & (tpm_est == 0)),
                TP_tpm = sum((tpm_truth > 0) * tpm_est),
                FP_tpm = sum((tpm_truth == 0) * tpm_est)) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
    
    print(nrow(exp_data_hap_exp))
    
    exp_data_exp_is_hap <- exp_data %>%
      filter(is_hap) %>% 
      select(-transcript, -HaplotypeProbability, -is_hap) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
    
    print(nrow(exp_data_exp_is_hap))
    
    print(exp_data_stats)
    
    dir_split <- strsplit(dirname(f), "/")[[1]]
    reads <- paste(dir_split[5], dir_split[6], sep = "_")
    method <- dir_split[7]
    graph <- dir_split[8]
    
    if (grepl("map_fast", basename(f))) {
      
        method = paste(method, "map", sep = "_")
    }

    exp_data_hap_prob <- exp_data_hap_prob %>%
      mutate(Reads = reads) %>%
      mutate(Method = method) %>%
      mutate(Graph = graph)  
    
    exp_data_hap_exp <- exp_data_hap_exp %>%
      mutate(Reads = reads) %>%
      mutate(Method = method) %>%
      mutate(Graph = graph)  
    
    exp_data_exp_is_hap <- exp_data_exp_is_hap %>%
      mutate(Reads = reads) %>%
      mutate(Method = method) %>%
      mutate(Graph = graph)  
    
    exp_data_stats <- exp_data_stats %>%
      mutate(Reads = reads) %>%
      mutate(Method = method) %>%
      mutate(Graph = graph)  
  
    save(exp_data_hap_prob, exp_data_hap_exp, exp_data_exp_is_hap, exp_data_stats, file = paste("rdata/", reads, method, graph, ".RData", sep = "", collapse = ""))
    
    gc()
  }
}

########

dataset <- "SRR1153470"
sim_mean <- 277
sim_sd <- 43

for (ref_name in c("gencode100", "1kg_NA12878_gencode100", "1kg_EURnonCEU_af002_gencode100", "1kg_nonCEU_af001_gencode100", "1kg_all_af001_gencode100")) {

  parse_data(dataset, sim_mean, sim_sd, "sim_vg", ref_name)
}

for (ref_name in c("1kg_EURnonCEU_af002_gencode100_unidi", "1kg_nonCEU_af001_gencode100_unidi", "1kg_all_af001_gencode100_unidi")) {

  parse_data(dataset, sim_mean, sim_sd, "real", ref_name)
}

dataset <- "ENCSR000AED_rep1"
sim_mean <- 216
sim_sd <- 24

for (ref_name in c("gencode100", "1kg_NA12878_gencode100", "1kg_EURnonCEU_af002_gencode100", "1kg_nonCEU_af001_gencode100", "1kg_all_af001_gencode100")) {

  parse_data(dataset, sim_mean, sim_sd, "sim_vg", ref_name)
}

for (ref_name in c("1kg_EURnonCEU_af002_gencode100_unidi", "1kg_nonCEU_af001_gencode100_unidi", "1kg_all_af001_gencode100_unidi")) {

  parse_data(dataset, sim_mean, sim_sd, "real", ref_name)
}

########
