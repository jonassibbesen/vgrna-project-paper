
# plot_hst_diversity.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("ggrepel")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


set.seed(1234)

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  data <- read_table2(filename, col_names = T)

  return(data)
}

wes_cols <- wes_palette("Cavalcanti1")

hst_data <- map_dfr(list.files(pattern=".*_hst_stats.txt", full.names = T, recursive = T), parse_file)

hst_data <- hst_data %>%
  group_by(name, pop, spop) %>%
  summarise(total = sum(total), num_sam = sum(num_sam), num_pop = sum(num_pop), num_spop = sum(num_spop)) %>%
  mutate("Excluding\nsample" = 1 - num_sam / total) %>%
  mutate("Excluding\npopulation" = 1 - num_pop / total) %>%
  mutate("Excluding\nsuper popluation" = 1 - num_spop / total) %>%
  gather("Excluding\nsample", "Excluding\npopulation", "Excluding\nsuper popluation", key = "type", value = "value") 

hst_data$type <- factor(hst_data$type, levels = c("Excluding\nsample", "Excluding\npopulation", "Excluding\nsuper popluation"))


pdf("plots/hst_diversity.pdf", height = 6, width = 9, pointsize = 12)
hst_data %>%
  ggplot(aes(x = type, y = value, fill = spop)) +
  geom_boxplot() +
  scale_fill_manual(values = wes_cols) +
  labs(fill = "Super\npopulation", x = "") +
  ylab("Fraction unique haplotype-specific transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16)) 
dev.off()

