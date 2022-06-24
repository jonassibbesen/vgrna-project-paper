
# plot_gibbs_stats.R

rm(list=ls())

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_r2/")

########

gibbs_samples <- read_table("gibbs/rpvg_gibbs_mpmap_1kg_nonCEU_af001_gencode100_sim_vg_r2_ENCSR000AED_rep1_gibbs_stats.txt") %>%
  select(-ClusterID)

exp_sim <- read_table("truth/truth_exp_sim_vg_ENCSR000AED_rep1_1kg_nonCEU_af001_gencode100.tsv") %>%
  select(-tpm_truth) %>%
  rename(Name = name)

gibbs_samples_sim <- gibbs_samples %>% 
  left_join(exp_sim, by = "Name") %>%
  replace_na(list(count_truth = 0 )) 

gibbs_samples_sim %>%
  mutate(inCI90 = (count_truth >= CILower90 & count_truth <= CIUpper90)) %>%
  mutate(belowCI90 = (count_truth < CILower90)) %>%
  mutate(aboveCI90 = (count_truth > CIUpper90)) %>%
  summarise(fracInCI90 = sum(inCI90)/n(), fracBelowCI90 = sum(belowCI90)/n(), fracAboveCI90 = sum(aboveCI90)/n()) %>%
  print()

gibbs_samples_sim_exp <- gibbs_samples_sim %>%
  add_column(Reads = "Simulated reads (vg)") %>%
  add_column(FacetRow = "")

wes_cols <- c(wes_palette("FantasticFox1"))

max_val <- max(c(log10(gibbs_samples_sim_exp$count_truth + 1), log10(gibbs_samples_sim_exp$CIUpper90 + 1)))

png("plots/gibbs/sim_r2_gibbs_ci_scatter_ENCSR000AED_rep1_lower.png", width = 4, height = 4, units = 'in', res = 600, pointsize = 12)
gibbs_samples_sim_exp %>%
  ggplot(aes(x = log10(CILower90 + 1), y = log10(count_truth + 1))) +
  geom_point(size = 0.25, alpha = 0.15, col = wes_cols[3]) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Lower credibility interval (log10 + 1)") +
  ylab("Simulated expression (log10 + 1)") +
  facet_grid(FacetRow ~ Reads) +
  xlim(c(0, max_val)) +
  ylim(c(0, max_val)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

png("plots/gibbs/sim_r2_gibbs_ci_scatter_ENCSR000AED_rep1_upper.png", width = 4, height = 4, units = 'in', res = 600, pointsize = 12)
gibbs_samples_sim_exp %>%
  ggplot(aes(x = log10(CIUpper90 + 1), y = log10(count_truth + 1))) +
  geom_point(size = 0.25, alpha = 0.15, col = wes_cols[4]) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Upper credibility interval (log10 + 1)") +
  ylab("Simulated expression (log10 + 1)") +
  facet_grid(FacetRow ~ Reads) +
  xlim(c(0, max_val)) +
  ylim(c(0, max_val)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 12))
dev.off()

########
