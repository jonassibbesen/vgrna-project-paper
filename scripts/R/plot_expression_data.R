
# plot_expression_data.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_r2/")

########

set.seed(4321)

exp_data_hap_prob_all <- list()
exp_data_hap_exp_all <- list()
exp_data_exp_is_hap_all <- list()
exp_data_stats_all <- list()

prepareData <- function(data) {
  
  data$Reads = recode_factor(data$Reads,
                              "sim_vg_SRR1153470" = "Simulated reads (vg)",
                              "real_SRR1153470" = "Real reads",
                              "sim_vg_ENCSR000AED_rep1" = "Simulated reads (vg)",
                              "real_ENCSR000AED_rep1" = "Real reads",
                              "sim_vg_r1_SRR1153470" = "Simulated reads (vg)",
                              "real_r1_SRR1153470" = "Real reads",
                              "sim_vg_r1_ENCSR000AED_rep1" = "Simulated reads (vg)",
                              "real_r1_ENCSR000AED_rep1" = "Real reads",
                              "sim_vg_r2_SRR1153470" = "Simulated reads (vg)",
                              "real_r2_SRR1153470" = "Real reads",
                              "sim_vg_r2_ENCSR000AED_rep1" = "Simulated reads (vg)",
                              "real_r2_ENCSR000AED_rep1" = "Real reads",
                              "sim_vg_debug_SRR1153470" = "Simulated reads (vg)"
                             )

  data$Method = recode_factor(data$Method,
                              "kallisto" = "Kallisto",
                              "kallisto_strand" = "Kallisto",
                              "kallisto_strand_bias" = "Kallisto (bias)",
                              "salmon" = "Salmon",
                              "salmon_bias" = "Salmon (bias)",
                              "salmon_em" = "Salmon (EM)",
                              "salmon_em_bias" = "Salmon (EM, bias)",
                              "rsem" = "RSEM (def)",
                              "rsem_k1k" = "RSEM",
                              "rsem_k2k" = "RSEM (k2000)",
                              "rsem_strand" = "RSEM (def)",
                              "rsem_strand_k1k" = "RSEM",
                              "rsem_strand_k2k" = "RSEM (k2000)",
                              "rpvg" = "rpvg",
                              "rpvg_trs" = "rpvg (reference)",
                              "rpvg_gam" = "rpvg (vg mpmap, single-path)",
                              "rpvg_map" = "rpvg (vg map, single-path)",
                              "rpvg_strand" = "rpvg",
                              "rpvg_strand_gam" = "rpvg (vg mpmap, single-path)",
                              "rpvg_strand_map" = "rpvg (vg map, single-path)",
                              "rpvg_paper" = "rpvg_paper"
                              )
  
  data$Graph = recode_factor(data$Graph,
                             "gencode100" = "Reference",
                             "gencode100_genes" = "Reference (genes)",
                             "1kg_NA12878_gencode100" = "Personal (NA12878)",
                             "1kg_NA12878_gencode100_genes" = "Personal (NA12878, genes)",
                             "1kg_EURnonCEU_af002_gencode100" = "Europe (excl. CEU)",
                             "1kg_EURnonCEU_af002_gencode100_genes" = "Europe (excl. CEU, genes)",
                             "1kg_EURnonCEU_af002_gencode100_unidi" = "Europe (excl. CEU)",
                             "1kg_EURnonCEU_af002_gencode100_unidi_mt_sca" = "Europe (excl. CEU)",
                             "1kg_EURnonCEU_af002_gencode100_unidi_decoy_mt_sca" = "Europe (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100" = "Whole (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100_genes" ="Whole (excl. CEU, genes)",
                             "1kg_nonCEU_af001_gencode100_unidi" = "Whole (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100_unidi_mt_sca" = "Whole (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100_unidi_decoy_mt_sca" = "Whole (excl. CEU)",
                             "1kg_all_af001_gencode100" = "Whole",
                             "1kg_all_af001_gencode100_genes" = "Whole (genes)",
                             "1kg_all_af001_gencode100_unidi" = "Whole",
                             "1kg_all_af001_gencode100_unidi_mt_sca" = "Whole",
                             "1kg_all_af001_gencode100_unidi_decoy_mt_sca" = "Whole"
)

  data$Reads <- factor(data$Reads, levels = c("Simulated reads (vg)", "Real reads"))
  data$Method <- factor(data$Method, levels = c("Kallisto", "Kallisto (bias)", "Salmon", "Salmon (bias)", "Salmon (EM)", "Salmon (EM, bias)", "RSEM (def)", "RSEM", "RSEM (k2000)", "rpvg", "rpvg (reference)", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)", "rpvg_paper"))
  data$Graph <- factor(data$Graph, levels = c("Reference", "Reference (genes)", "Personal (NA12878)", "Personal (NA12878, genes)", "Europe (excl. CEU)", "Europe (excl. CEU, genes)", "Whole (excl. CEU)", "Whole (excl. CEU, genes)", "Whole", "Whole (genes)"))
  
  data <- data %>%
    rename(Pantranscriptome = Graph)

  return(data)
}

args = commandArgs(trailingOnly=TRUE)
print(args)

dataset <- args[1]

for (f in list.files(path = "./rdata", pattern = paste(".*", dataset, ".*RData", sep = ""), full.names = T, recursive = T)) { 
  
  print(f)
  load(f)
  
  exp_data_hap_prob_all[[f]] <- prepareData(exp_data_hap_prob)
  exp_data_hap_exp_all[[f]] <- prepareData(exp_data_hap_exp)
  exp_data_exp_is_hap_all[[f]] <- prepareData(exp_data_exp_is_hap)
  exp_data_stats_all[[f]] <- prepareData(exp_data_stats)
}

# for (f in list.files(path = "../quant_r1/rdata", pattern = paste(".*", dataset, ".*RData", sep = ""), full.names = T, recursive = T)) {
# 
#   if (!grepl("rpvg", f) | grepl("rpvg_amq", f) | grepl("rpvg_gam", f) | grepl("rpvg_map", f) | grepl("rpvg_strand_amq", f) | grepl("rpvg_strand_gam", f) | grepl("rpvg_strand_map", f)) {
# 
#     next
#   }
# 
#   print(f)
#   load(f)
# 
#   exp_data_hap_prob$Method = recode_factor(exp_data_hap_prob$Method,
#                                            "rpvg" = "rpvg_paper",
#                                            "rpvg_strand" = "rpvg_paper")
# 
#   exp_data_hap_exp$Method = recode_factor(exp_data_hap_exp$Method,
#                                           "rpvg" = "rpvg_paper",
#                                           "rpvg_strand" = "rpvg_paper")
# 
#   exp_data_exp_is_hap$Method = recode_factor(exp_data_exp_is_hap$Method,
#                                           "rpvg" = "rpvg_paper",
#                                           "rpvg_strand" = "rpvg_paper")
# 
#   exp_data_stats$Method = recode_factor(exp_data_stats$Method,
#                                         "rpvg" = "rpvg_paper",
#                                         "rpvg_strand" = "rpvg_paper")
# 
#   exp_data_hap_prob_all[[f]] <- prepareData(exp_data_hap_prob)
#   exp_data_hap_exp_all[[f]] <- prepareData(exp_data_hap_exp)
#   exp_data_exp_is_hap_all[[f]] <- prepareData(exp_data_exp_is_hap)
#   exp_data_stats_all[[f]] <- prepareData(exp_data_stats)
# }

do.call(bind_rows, exp_data_hap_exp_all) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  filter(tpm_est > 0) %>%
  summarise(min_tpm_est = min(tpm_est), max_tpm_est = max(tpm_est), sum_exp = sum(tpm_est > 0)) %>%
  print(n = 100)

########

exp_data_hap_prob_all_roc <- do.call(bind_rows, exp_data_hap_prob_all) %>%
  mutate(HaplotypeProbability = signif(HaplotypeProbability, 2)) %>%
  group_by(HaplotypeProbability, Reads, Method, Pantranscriptome) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  arrange(desc(HaplotypeProbability), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_tpm = TP_cs / max(TP_cs + FN_cs)) %>%
  mutate(PPV_tpm = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) 
  
exp_data_hap_prob_all_roc_points <- exp_data_hap_prob_all_roc %>% 
  mutate(HaplotypeProbability_floor = floor(HaplotypeProbability * 20) / 20) %>%
  group_by(HaplotypeProbability_floor, Reads, Method, Pantranscriptome) %>%
  mutate(is_min = (HaplotypeProbability == min(HaplotypeProbability))) %>%
  filter(is_min) %>%
  mutate(HaplotypeProbability = HaplotypeProbability_floor) %>%
  filter(HaplotypeProbability_floor == 0 | HaplotypeProbability_floor == 0.8 | HaplotypeProbability_floor == 1)


exp_data_hap_exp_all_roc <- do.call(bind_rows, exp_data_hap_exp_all) %>%
  mutate(tpm_est_sig = signif(tpm_est, 1)) %>%
  mutate(tpm_est_sig = ifelse(tpm_est_sig >= 10, 10, tpm_est_sig)) %>%
  group_by(tpm_est_sig, Reads, Method, Pantranscriptome) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  arrange(desc(tpm_est_sig), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_tpm = TP_cs / max(TP_cs + FN_cs)) %>%
  mutate(PPV_tpm = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) %>%
  filter(tpm_est_sig > 0)

exp_data_hap_exp_all_roc_points <- exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  mutate(is_min = (tpm_est_sig == min(tpm_est_sig))) %>%
  mutate(is_max = (tpm_est_sig == max(tpm_est_sig))) %>%
  filter(tpm_est_sig == 0.1 | tpm_est_sig == 1 | is_max) 


exp_data_stats_all_scatter <- do.call(bind_rows, exp_data_exp_is_hap_all) %>%
  filter(Reads == "Simulated reads (vg)") 


exp_data_stats_all_bars <- do.call(bind_rows, exp_data_stats_all) 

########

wes_cols <- c(wes_palette("GrandBudapest1")[1], wes_palette("GrandBudapest2")[4], wes_palette("GrandBudapest1")[2], wes_palette("Chevalier1")[c(1,2)], wes_palette("Darjeeling2")[c(2,3,4,5)])

exp_data_hap_prob_all_roc$FacetRow <- ""
exp_data_hap_prob_all_roc_points$FacetRow <- ""

exp_data_hap_prob_all_roc_points_sim <- exp_data_hap_prob_all_roc_points %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/prob/sim_r2_hap_prob_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_prob_all_roc %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = HaplotypeProbability)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_prob_all_roc_points_sim, size = 2) +
  geom_label_repel(data = exp_data_hap_prob_all_roc_points_sim[exp_data_hap_prob_all_roc_points_sim$Pantranscriptome == "Whole (excl. CEU)",], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
  scale_color_manual(values = wes_cols[c(4)]) +
  facet_grid(FacetRow ~ Reads) +
  xlab("Transcript expression precision") +
  ylab("Transcript expression recall") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 13))
dev.off()

exp_data_hap_prob_all_roc_points_real <- exp_data_hap_prob_all_roc_points %>%
  filter(Reads == "Real reads") %>%
  filter(Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/prob/real_r2_hap_prob_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_prob_all_roc %>%
  filter(Reads == "Real reads") %>%
  filter(Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = HaplotypeProbability)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_prob_all_roc_points_real, size = 2) +
  geom_label_repel(data = exp_data_hap_prob_all_roc_points_real[exp_data_hap_prob_all_roc_points_real$Pantranscriptome == "Whole (excl. CEU)",], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
  scale_color_manual(values = wes_cols[c(4)]) +
  facet_grid(FacetRow ~ Reads) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("Number of expressed non-NA12878 transcripts") +
  ylab("Number of expressed NA12878 transcripts") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 13))
dev.off()


exp_data_hap_exp_all_roc$FacetRow <- ""
exp_data_hap_exp_all_roc_points$FacetRow <- ""

exp_data_hap_exp_all_roc_points_sim <- exp_data_hap_exp_all_roc_points %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/hap/sim_r2_expression_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_exp_all_roc %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_sim, size = 2) +
  geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim[exp_data_hap_exp_all_roc_points_sim$Pantranscriptome == "Whole (excl. CEU)",], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
  scale_color_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("Transcript expression precision") +
  ylab("Transcript expression recall") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 14))
dev.off()

exp_data_hap_exp_all_roc_points_sim_sample <- exp_data_hap_exp_all_roc_points %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/hap/sim_r2_expression_roc_pre_sample_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_exp_all_roc %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_sim_sample, size = 2) +
  geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim_sample, size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
  scale_color_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("Transcript expression precision") +
  ylab("Transcript expression recall") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 14))
dev.off()

exp_data_hap_exp_all_roc_points_real <- exp_data_hap_exp_all_roc_points %>%
  filter(Reads == "Real reads") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/hap/real_r2_expression_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_exp_all_roc %>%
  filter(Reads == "Real reads") %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig, segment.color = Method)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_real, size = 2) +
  geom_label_repel(data = exp_data_hap_exp_all_roc_points_real[exp_data_hap_exp_all_roc_points_real$Pantranscriptome == "Whole (excl. CEU)",], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
  scale_color_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("Number of expressed non-NA12878 transcripts") +
  ylab("Number of expressed NA12878 transcripts") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 13))
dev.off()


exp_data_hap_exp_all_roc_rpvg <- exp_data_hap_exp_all_roc %>%
  filter(Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)")

exp_data_hap_exp_all_roc_points_rpvg <- exp_data_hap_exp_all_roc_points %>%
  filter(Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)")

exp_data_hap_exp_all_roc_rpvg$Method = recode_factor(exp_data_hap_exp_all_roc_rpvg$Method, "rpvg" = "rpvg (vg mpmap)")
exp_data_hap_exp_all_roc_points_rpvg$Method = recode_factor(exp_data_hap_exp_all_roc_points_rpvg$Method, "rpvg" = "rpvg (vg mpmap)")

exp_data_hap_exp_all_roc_rpvg$Method <- factor(exp_data_hap_exp_all_roc_rpvg$Method, levels = c("rpvg (vg mpmap)", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))
exp_data_hap_exp_all_roc_points_rpvg$Method <- factor(exp_data_hap_exp_all_roc_points_rpvg$Method, levels = c("rpvg (vg mpmap)", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

exp_data_hap_exp_all_roc_points_rpvg_sim <- exp_data_hap_exp_all_roc_points_rpvg %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/hap/sim_r2_expression_roc_pre_rpvg_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_exp_all_roc_rpvg %>%
  filter(Reads == "Simulated reads (vg)") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_rpvg_sim, size = 2) +
  geom_label_repel(data = exp_data_hap_exp_all_roc_points_rpvg_sim[exp_data_hap_exp_all_roc_points_rpvg_sim$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_rpvg_sim$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
  scale_color_manual(values = wes_cols[c(4,5,6)]) +
  facet_grid(FacetRow ~ Reads) +
  xlab("Transcript expression precision") +
  ylab("Transcript expression recall") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 13))
dev.off()

exp_data_hap_exp_all_roc_points_rpvg_real <- exp_data_hap_exp_all_roc_points_rpvg %>%
  filter(Reads == "Real reads") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole")

pdf(paste("plots/hap/real_r2_expression_roc_pre_rpvg_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_hap_exp_all_roc_rpvg %>%
  filter(Reads == "Real reads") %>%
  filter(Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_rpvg_real, size = 2) +
  geom_label_repel(data = exp_data_hap_exp_all_roc_points_rpvg_real[exp_data_hap_exp_all_roc_points_rpvg_real$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_rpvg_real$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
  scale_color_manual(values = wes_cols[c(4,5,6)]) +
  facet_grid(FacetRow ~ Reads) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("Number of expressed non-NA12878 transcripts") +
  ylab("Number of expressed NA12878 transcripts") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 12))
dev.off()

  
exp_data_stats_all_scatter$Method <- factor(exp_data_stats_all_scatter$Method, levels = c("rpvg", "Kallisto", "Salmon", "RSEM", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

pdf(paste("plots/exp/sim_r2_exp_scatter_hap_tpm_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
exp_data_stats_all_scatter %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Personal (NA12878)" | Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  mutate(value_x_log = log10(tpm_est + 1)) %>%
  mutate(value_y_log = log10(tpm_truth + 1)) %>%
  ggplot(aes(y = value_y_log, x = value_x_log)) +
  geom_bin2d(binwidth = c(0.1, 0.1)) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 0.25, alpha = 0.5) +
  facet_grid(Pantranscriptome ~ Method) +
  scale_fill_gradient(name = "Count", trans = "log10") +
  xlab("Estimated expression (log10 + 1)") +
  ylab("Simulated expression (log10 + 1)") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.2, "cm")) +
  theme(legend.key.width = unit(1.2, "cm")) +
  theme(text = element_text(size = 8)) +
  theme(legend.text = element_text(size = 8))
dev.off()


exp_data_stats_all_bars_main <- exp_data_stats_all_bars %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)") %>%
  filter(Pantranscriptome == "Reference" | Pantranscriptome == "Personal (NA12878)" | Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  ungroup() %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "All", frac_hap_error_tpm = 1) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", frac_hap_error_tpm = 1) %>%
  add_row(Reads = "Real reads", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "All", frac_hap_error_tpm = 1) %>%
  add_row(Reads = "Real reads", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", frac_hap_error_tpm = 1)

exp_data_stats_all_bars_main$Reads <- factor(exp_data_stats_all_bars_main$Reads, levels = c("Simulated reads (vg)", "Real reads"))
exp_data_stats_all_bars_main$Method <- factor(exp_data_stats_all_bars_main$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

exp_data_stats_all_bars_main$Pantranscriptome = recode_factor(exp_data_stats_all_bars_main$Pantranscriptome,
                           "Reference" = "Reference",
                           "Personal (NA12878)" = "Personal\n(NA12878)",
                           "Europe (excl. CEU)" = "Europe\n(excl. CEU)",
                           "Whole (excl. CEU)" = "Whole\n(excl. CEU)",
                           "Whole" = "Whole")

exp_data_stats_all_bars_main$Pantranscriptome <- factor(exp_data_stats_all_bars_main$Pantranscriptome, levels = c("Reference", "Personal\n(NA12878)", "Europe\n(excl. CEU)", "Whole\n(excl. CEU)", "Whole"))

exp_data_stats_all_bars_main$FacetRow <- ""

pdf(paste("plots/hap/sim_r2_real_hap_tpm_error_", dataset, ".pdf", sep = ""), height = 4, width = 7, pointsize = 12)
exp_data_stats_all_bars_main %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Europe\n(excl. CEU)" | Pantranscriptome == "Whole\n(excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = 1 - frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Fraction TPM on NA12878 haplotypes") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 11))
dev.off()


exp_data_stats_all_bars_main_rpvg <- exp_data_stats_all_bars_main %>%
  filter(Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)")

exp_data_stats_all_bars_main_rpvg$Method = recode_factor(exp_data_stats_all_bars_main_rpvg$Method, "rpvg" = "rpvg (vg mpmap)")
exp_data_stats_all_bars_main_rpvg$Method <- factor(exp_data_stats_all_bars_main_rpvg$Method, levels = c("rpvg (vg mpmap)", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

pdf(paste("plots/hap/sim_r2_real_hap_tpm_error_rpvg_", dataset, ".pdf", sep = ""), height = 4, width = 7.5, pointsize = 12)
exp_data_stats_all_bars_main_rpvg %>%
  filter(Pantranscriptome == "Europe\n(excl. CEU)" | Pantranscriptome == "Whole\n(excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = 1 - frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols[c(4,5,6)]) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Fraction TPM on NA12878 haplotypes") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 11))
dev.off()


exp_data_stats_all_bars_main_sim <- exp_data_stats_all_bars %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "RSEM (def)" | Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)") %>%
  filter(Pantranscriptome == "Reference" | Pantranscriptome == "Personal (NA12878)" | Pantranscriptome == "Europe (excl. CEU)" | Pantranscriptome == "Whole (excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ungroup() %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "All", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "Expressed", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "Haplotype", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole (excl. CEU)", Type = "Transcript", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole", Type = "Expressed", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole", Type = "Haplotype", Spearman_tpm = 0, ARD_mean_tpm = 0) %>%
  add_row(Reads = "Simulated reads (vg)", Method = "RSEM", Pantranscriptome = "Whole", Type = "Transcript", Spearman_tpm = 0, ARD_mean_tpm = 0)

exp_data_stats_all_bars_main_sim$Type <- as.factor(exp_data_stats_all_bars_main_sim$Type)
exp_data_stats_all_bars_main_sim$Method <- factor(exp_data_stats_all_bars_main_sim$Method, levels = c("Kallisto", "Salmon", "RSEM", "RSEM (def)", "rpvg", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

exp_data_stats_all_bars_main_sim$Pantranscriptome = recode_factor(exp_data_stats_all_bars_main_sim$Pantranscriptome,
                                                              "Reference" = "Reference",
                                                              "Personal (NA12878)" = "Personal\n(NA12878)",
                                                              "Europe (excl. CEU)" = "Europe\n(excl. CEU)",
                                                              "Whole (excl. CEU)" = "Whole\n(excl. CEU)",
                                                              "Whole" = "Whole")

exp_data_stats_all_bars_main_sim$Pantranscriptome <- factor(exp_data_stats_all_bars_main_sim$Pantranscriptome, levels = c("Reference", "Personal\n(NA12878)", "Europe\n(excl. CEU)", "Whole\n(excl. CEU)", "Whole"))

exp_data_stats_all_bars_main_sim$FacetRow <- ""


exp_data_stats_all_bars_main_sim_all_hap <- exp_data_stats_all_bars_main_sim %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Personal\n(NA12878)" | Pantranscriptome == "Europe\n(excl. CEU)" | Pantranscriptome == "Whole\n(excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Type == "All" | Type == "Haplotype")

exp_data_stats_all_bars_main_sim_all_hap$Type <- factor(exp_data_stats_all_bars_main_sim_all_hap$Type, levels = c("All", "Haplotype"))

pdf(paste("plots/exp/sim_r2_mard_all_hap_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_all_hap, Type == "All"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_all_hap, Type == "Haplotype"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Mean relative expression difference") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()

pdf(paste("plots/exp/sim_r2_spearman_all_hap_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_all_hap, Type == "All"), aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_all_hap, Type == "Haplotype"), aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Spearman expression correlation") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()


exp_data_stats_all_bars_main_sim_rpvg <- exp_data_stats_all_bars_main_sim %>%
  filter(Method == "rpvg" | Method == "rpvg (vg mpmap, single-path)" | Method == "rpvg (vg map, single-path)") %>%
  filter(Pantranscriptome == "Personal\n(NA12878)" | Pantranscriptome == "Europe\n(excl. CEU)" | Pantranscriptome == "Whole\n(excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Type == "All" | Type == "Haplotype")

exp_data_stats_all_bars_main_sim_rpvg$Type <- factor(exp_data_stats_all_bars_main_sim_rpvg$Type, levels = c("All", "Haplotype"))

exp_data_stats_all_bars_main_sim_rpvg$Method = recode_factor(exp_data_stats_all_bars_main_sim_rpvg$Method, "rpvg" = "rpvg (vg mpmap)")
exp_data_stats_all_bars_main_sim_rpvg$Method <- factor(exp_data_stats_all_bars_main_sim_rpvg$Method, levels = c("rpvg (vg mpmap)", "rpvg (vg mpmap, single-path)", "rpvg (vg map, single-path)"))

pdf(paste("plots/exp/sim_r2_mard_all_hap_tpm_rpvg_", dataset, ".pdf", sep = ""), height = 4.5, width = 7.5, pointsize = 12)
ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_rpvg, Type == "All"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_main_sim_rpvg, Type == "Haplotype"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols[c(4,5,6)]) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Mean relative expression difference") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12.5))
dev.off()


exp_data_stats_all_bars_main_sim_expressed <- exp_data_stats_all_bars_main_sim %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
  filter(Pantranscriptome == "Personal\n(NA12878)" | Pantranscriptome == "Europe\n(excl. CEU)" | Pantranscriptome == "Whole\n(excl. CEU)" | Pantranscriptome == "Whole") %>%
  filter(Type == "Expressed")

exp_data_stats_all_bars_main_sim_expressed$Type <- factor(exp_data_stats_all_bars_main_sim_expressed$Type, levels = c("Expressed"))

pdf(paste("plots/exp/sim_r2_mard_exp_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
ggplot() +
  geom_bar(data = exp_data_stats_all_bars_main_sim_expressed, aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  scale_alpha_manual(name = "Transcripts", values = c(1), labels = c("Expressed"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Mean relative expression difference") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()

pdf(paste("plots/exp/sim_r2_spearman_exp_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
ggplot() +
  geom_bar(data = exp_data_stats_all_bars_main_sim_expressed, aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  scale_alpha_manual(name = "Transcripts", values = c(1), labels = c("Expressed"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Spearman expression correlation") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()


exp_data_stats_all_bars_main_sim_transcript <- exp_data_stats_all_bars_main_sim %>%
  filter(Method == "Kallisto" | Method == "Salmon" | (Method == "RSEM (def)" & Pantranscriptome == "Reference") | (Method == "RSEM" & Pantranscriptome != "Reference") | Method == "rpvg") %>%
  filter(Type == "Transcript")

exp_data_stats_all_bars_main_sim_transcript$Method = recode_factor(exp_data_stats_all_bars_main_sim_transcript$Method, "RSEM (def)" = "RSEM")
exp_data_stats_all_bars_main_sim_transcript$Method <- factor(exp_data_stats_all_bars_main_sim_transcript$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg"))

exp_data_stats_all_bars_main_sim_transcript$Type <- factor(exp_data_stats_all_bars_main_sim_transcript$Type, levels = c("Transcript"))

pdf(paste("plots/exp/sim_r2_mard_ref_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 7, pointsize = 12)
ggplot() +
  geom_bar(data = exp_data_stats_all_bars_main_sim_transcript, aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  scale_alpha_manual(name = "Transcripts", values = c(1), labels = c("Reference"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Mean relative expression difference") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()

pdf(paste("plots/exp/sim_r2_spearman_ref_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 7, pointsize = 12)
ggplot() +
  geom_bar(data = exp_data_stats_all_bars_main_sim_transcript, aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  scale_alpha_manual(name = "Transcripts", values = c(1), labels = c("Reference"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ Reads) +
  xlab("") +
  ylab("Spearman expression correlation") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 12))
dev.off()

######## Debug

print("DEBUG")

exp_data_stats_all_bars %>%
  filter(Pantranscriptome != "NA12878") %>%
  filter(Type == "All") %>%
  select(Reads, Method, Pantranscriptome, n, frac_hap_error_count, frac_hap_error_tpm) %>%
  print(n = 100)

exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  filter(tpm_est_sig == 0.01) %>%
  summarise(max_FP_cs = max(FP_cs)) %>%
  print(n = 100)


pdf(paste("plots/debug/hap_prob_rocs_pt_debug_", dataset, ".pdf", sep = ""), width = 9, height = 9)

for (pt in unique(exp_data_hap_prob_all_roc$Pantranscriptome)) {

  exp_data_hap_prob_all_roc_points_pt <- exp_data_hap_prob_all_roc_points %>%
    filter(Pantranscriptome == pt)

  delta <- max(max(exp_data_hap_prob_all_roc$PPV_tpm) - min(exp_data_hap_prob_all_roc$PPV_tpm), max(exp_data_hap_prob_all_roc$TPR_tpm) - min(exp_data_hap_prob_all_roc$TPR_tpm))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$PPV_tpm), min(exp_data_hap_prob_all_roc$PPV_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$TPR_tpm), min(exp_data_hap_prob_all_roc$TPR_tpm) + delta)) +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression recall (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  delta <- max(max(exp_data_hap_prob_all_roc$FP_cs) - min(exp_data_hap_prob_all_roc$FP_cs), max(exp_data_hap_prob_all_roc$TP_cs) - min(exp_data_hap_prob_all_roc$TP_cs))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    scale_x_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$FP_cs), min(exp_data_hap_prob_all_roc$FP_cs) + delta)) +
    scale_y_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$TP_cs), min(exp_data_hap_prob_all_roc$TP_cs) + delta)) +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  delta <- max(max(exp_data_hap_prob_all_roc$frac_incorrect_tpm) - min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), max(exp_data_hap_prob_all_roc$frac_correct_tpm) - min(exp_data_hap_prob_all_roc$frac_correct_tpm))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), min(exp_data_hap_prob_all_roc$frac_incorrect_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$frac_correct_tpm), min(exp_data_hap_prob_all_roc$frac_correct_tpm) + delta)) +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/debug/hap_prob_rocs_me_debug_", dataset, ".pdf", sep = ""), width = 9, height = 9)

for (me in unique(exp_data_hap_prob_all_roc$Method)) {

  exp_data_hap_prob_all_roc_points_me <- exp_data_hap_prob_all_roc_points %>%
    filter(Method == me)

  delta <- max(max(exp_data_hap_prob_all_roc$PPV_tpm) - min(exp_data_hap_prob_all_roc$PPV_tpm), max(exp_data_hap_prob_all_roc$TPR_tpm) - min(exp_data_hap_prob_all_roc$TPR_tpm))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$PPV_tpm), min(exp_data_hap_prob_all_roc$PPV_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$TPR_tpm), min(exp_data_hap_prob_all_roc$TPR_tpm) + delta)) +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression recall (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  delta <- max(max(exp_data_hap_prob_all_roc$FP_cs) - min(exp_data_hap_prob_all_roc$FP_cs), max(exp_data_hap_prob_all_roc$TP_cs) - min(exp_data_hap_prob_all_roc$TP_cs))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    scale_x_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$FP_cs), min(exp_data_hap_prob_all_roc$FP_cs) + delta)) +
    scale_y_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$TP_cs), min(exp_data_hap_prob_all_roc$TP_cs) + delta)) +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  delta <- max(max(exp_data_hap_prob_all_roc$frac_incorrect_tpm) - min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), max(exp_data_hap_prob_all_roc$frac_correct_tpm) - min(exp_data_hap_prob_all_roc$frac_correct_tpm))

  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), min(exp_data_hap_prob_all_roc$frac_incorrect_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$frac_correct_tpm), min(exp_data_hap_prob_all_roc$frac_correct_tpm) + delta)) +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/debug/expression_rocs_pt_debug_", dataset, ".pdf", sep = ""))

for (pt in unique(exp_data_hap_exp_all_roc$Pantranscriptome)) {

  exp_data_hap_exp_all_roc_points_pt <- exp_data_hap_exp_all_roc_points %>%
    filter(Pantranscriptome == pt)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression recall (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/debug/expression_rocs_me_debug_", dataset, ".pdf", sep = ""))

for (me in unique(exp_data_hap_exp_all_roc$Method)) {

  exp_data_hap_exp_all_roc_points_me <- exp_data_hap_exp_all_roc_points %>%
    filter(Method == me)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression recall (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)

  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/debug/general_stats_debug_", dataset, ".pdf", sep = ""))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = frac_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Fraction estimated CPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = num_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = num_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = ExpCorrect, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Fraction correct expressed") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = Pearson_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = Pearson_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = Spearman_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Mean absolute difference (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads (vg)") %>%
  ggplot(aes(x = Pantranscriptome, y = ARD_mean_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Mean absolute difference (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

dev.off()

########
