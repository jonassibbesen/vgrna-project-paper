
# utils.R

library("dplyr")
library("readr")
library("tibble")
library("tidyr")
library("purrr")
library("stringr")
library("ggrepel")
library("scales")
library("truncnorm")
library("ggplot2")
library("wesanderson")

########

set.seed(4321)
options(scipen = 999)

wes_cols <- c(rev(wes_palette("Rushmore1")[c(1,3,4,5)]), wes_palette("Darjeeling1")[c(5)], wes_palette("Royal2")[c(3,2,1)])

########

printHeader <- function() {

  args <- commandArgs()
  cript_dir <- dirname(sub("--file=", "", args[4]))
  print(script_dir)
  
  print(args)
  system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
  system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
}

plotRocCurveMapq <- function(data, cols, plot_numbers, lt_title) {
  
  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
  data_roc <- data %>% 
    mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
    group_by(Method, LineType, FacetRow, FacetCol, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    group_by(Method, LineType, FacetRow, FacetCol) %>%
    mutate(N = sum(TP) + sum(FP)) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    group_by(LineType, FacetRow, FacetCol) %>%
    mutate(FNcs = N - TPcs - FPcs) %>%
    mutate(TPR = TPcs / (TPcs + FNcs), FDR = FPcs / (TPcs + FPcs)) %>%
    filter(MapQ >= 0)
  
  min_lim_x <- min(data_roc$TPR)
  
  data_roc %>% filter(MapQ == 60 | MapQ == 255) %>% print(n = 100)
  
  a <- annotation_logticks(sides = "l")
  a$data <- data.frame(x = NA, FacetCol = c(as.character(data_roc$FacetCol[1])))
  
  p <- data_roc %>%
    ggplot(aes(y = FDR, x = TPR, color = Method, linetype = LineType, shape = LineType, label = MapQ)) +
    a +
    geom_line(size = 1) +
    geom_point(data = subset(data_roc, MapQ == 0 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60 | MapQ == 255), size = 2, alpha = 1)
    
  if (plot_numbers) {
    
    p <- p +
      geom_text_repel(data = subset(data_roc, MapQ == 0 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60| MapQ == 255), size = 3, fontface = 2, show.legend = FALSE)
  }

  p <- p +
    scale_y_continuous(trans = 'log10') +
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    scale_linetype_discrete(name = lt_title) +
    scale_shape_discrete(name = lt_title) +
    xlim(c(min_lim_x, 1)) +
    xlab("Mapping recall") +
    ylab("Mapping error (1 - precision)") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 13)) +
    theme(legend.text = element_text(size = 12)) 
  print(p)   
}

plotMapQCurve <- function(data, cols, ylab) {

  p <- data %>%
    ggplot(aes(y = Value, x = MapQ, color = Method, linetype = Reference, shape = Reference)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    xlab("Mapping quality threshold") +
    ylab(ylab) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 13)) +
    theme(legend.text = element_text(size = 12))
  print(p)
}

plotErrorCurve <- function(data, cols) {
  
  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct) %>%
    filter(MapQ > 0) %>%
    mutate(MapQ = ifelse(MapQ > 60, 60, MapQ)) %>%
    group_by(Method, Reference, Reads, FacetRow, FacetCol, MapQ) %>%
    summarise(Count = sum(Count), TP = sum(TP), FP = sum(FP)) %>% 
    mutate(Est_MapQ = -10 * log10(FP / (TP +FP)))
  
  p <- data %>%
    ggplot(aes(y = Est_MapQ, x = MapQ, color = Method, size = Count)) +
    geom_point() +
    geom_line(size = 1) +
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    scale_size(breaks = c(0.5*10^6, 10*10^6, 75*10^6), labels = comma) +
    xlab("Estimated mapping quality") +
    ylab("Empirical mapping quality") +
    xlim(c(0,60)) +
    ylim(c(0,60)) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 13)) +
    theme(legend.text = element_text(size = 12))
  print(p) 
}

plotLogBins <- function(data, cols, xlab, ylab) {
  
  data <- data %>%
    mutate(value_x_log = log10(value_x + 1)) %>%
    mutate(value_y_log = log10(value_y + 1))

  max_val <- max(data$value_x_log, data$value_y_log)
  breaks <- c(1, 10, 100, 1000, 10000)
    
  p <- data %>%
    ggplot(aes(y = value_y_log, x = value_x_log)) +
    geom_bin2d() +
    facet_grid(FacetRow ~ FacetCol) +
    scale_fill_gradient(name = "Count", trans = "log10", breaks = breaks, labels = breaks) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(0, max_val)) +
    ylim(c(0, max_val)) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 10)) +
    theme(legend.text = element_text(size = 9))
  print(p) 
}

plotBiasCurve <- function(data, cols, min_count) {
  
  data <- data %>%
    mutate(ref = (ref_up + ref_down) / 2) %>%
    mutate(alt = (alt_up + alt_down) / 2) %>%
    filter(ref + alt >= min_count) %>%
    mutate(frac = alt / (ref + alt)) %>%
    mutate(len = ifelse(len > 15, 16, len)) %>%
    mutate(len = ifelse(len < -15, -16, len)) %>%
    group_by(Method, Reference, Reads, FacetCol, FacetRow, var, len) %>%
    summarise(n = n(), ref_count = sum(ref), alt_count = sum(alt), frac_mean = mean(frac))  
  
  p <- data %>% 
    ggplot(aes(y = frac_mean, x = len, color = Method, linetype = Reference, shape = Reference, label = sprintf("%0.3f", round(frac_mean, digits = 3)))) +
    geom_hline(yintercept = 0.5, size = 0.5, linetype = 1, alpha = 0.75) + 
    geom_line(size = 0.75) + 
    geom_point(data = subset(data, len == 0), size = 2) +  
    geom_text_repel(data = subset(data, len == 0), size = 3, fontface = 2, box.padding = 0.75, show.legend = FALSE) +  
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    scale_x_continuous(breaks=c(-16, -10, -5, 0, 5, 10, 16), labels = c("<-15", "-10", "-5", "SNV", "5", "10", ">15")) +
    ylim(c(0.3, 0.7)) +
    xlab("Allele length") +
    ylab("Mean fraction of reads on alt allele") +
    guides(linetype = "none") +
    guides(shape = "none") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.3, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 14)) +
    theme(legend.text = element_text(size = 13))
  print(p)
}

plotBiasBinom <- function(data, cols, alpha_thres, min_count) {
  
  data <- data %>%
    ungroup() %>%
    mutate(ref = (ref_up + ref_down) / 2) %>%
    mutate(alt = (alt_up + alt_down) / 2) %>%
    filter(ref + alt >= min_count) %>%
    mutate(count = ref + alt) %>%
    rowwise() %>%
    mutate(p_val = binom.test(x = c(round(ref), round(alt)), alternative = "two.sided")$p.value) %>%
    group_by(Method, Reference, Reads, FacetCol, FacetRow) %>%
    summarise(n = n(), n_bias = sum(p_val <= alpha_thres))
  
  data %>% print(n = 100)
  max_y <- max(data$n_bias / data$n)
  
  p <- data %>% 
    ggplot(aes(y = n_bias / n, x = n, color = Method, shape = Reference)) +
    geom_hline(yintercept = alpha_thres, size = 0.25, linetype = 1, alpha = 0.75) + 
    geom_point(size = 1.5) +
    scale_color_manual(values = cols) +
    facet_grid(FacetRow ~ FacetCol, scales="free_x") +
    ylim(c(0, max_y)) +
    xlab(bquote("Number of variants (coverage">=.(min_count)*")")) +
    ylab(bquote("Fraction of biased variants ("*alpha==.(alpha_thres)*")")) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(legend.key.width = unit(1.2, "cm")) +
    theme(text = element_text(size = 10)) +
    theme(legend.text = element_text(size = 9))
  print(p)
}

plotStatsBarTwoLayer <- function(data, cols) {
  
  data_bar <- data
  min_frac <- data_bar %>% ungroup() %>% filter(Frac > 0) %>% summarise(frac = min(Frac))
  
  p <- ggplot() +
    geom_bar(data = data_bar[data_bar$Filter == "MapQ >= 30",], aes(Reference, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge()) +
    geom_bar(data = data_bar[data_bar$Filter == "Mapped",], aes(Reference, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
    scale_fill_manual(values = cols) +
    scale_alpha_manual(name = "Filter", values = c(0.5, 1), labels = c("Mapped", bquote("MapQ">="30")), drop = F) + 
    scale_y_continuous(limits = c(floor(min_frac$frac * 10) / 10, 1), oob = rescale_none) +
    facet_grid(FacetRow ~ FacetCol) +
    xlab("") +
    ylab("Mapping rate") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(text = element_text(size = 12)) +
    theme(legend.text = element_text(size = 11))
  print(p)
}

plotStatsBar <- function(data, cols, ylab) {
  
  min_frac <- data %>% filter(Frac > 0) %>% ungroup() %>% summarise(frac = min(Frac))
  
  p <- data %>% 
    ggplot(aes(Reference, y = Frac, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(floor(min_frac$frac * 2) / 2, 1), oob = rescale_none) +
    facet_grid(FacetRow ~ FacetCol) +
    xlab("") +
    ylab(ylab) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(text = element_text(size = 12)) +
    theme(legend.text = element_text(size = 11))
  print(p)
}

plotBar <- function(data, cols, ylab) {
  
  p <- data %>% 
    ggplot(aes(x = Reference, y = Value, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    scale_fill_manual(values = cols) +
    facet_grid(FacetRow ~ FacetCol) +
    xlab("") +
    ylab(ylab) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(panel.spacing = unit(0.5, "cm")) +
    theme(text = element_text(size = 12)) +
    theme(legend.text = element_text(size = 11))
  print(p)
}

plotRocBenchmarkMapQ <- function(data, cols, lt_title, filename) {

  pdf(paste(filename, "_roc.pdf", sep = ""), height = 5, width = 7, pointsize = 12)
  plotRocCurveMapq(data, cols, T, lt_title)
  dev.off() 
}

plotRocBenchmarkMapQDebug <- function(data, cols, lt_title, filename) {
  
  pdf(paste(filename, "_roc.pdf", sep = ""), height = 7, width = 9, pointsize = 12)
  plotRocCurveMapq(data, cols, F, lt_title)
  dev.off() 
}

plotErrorBenchmark <- function(data, cols, filename) {
  
  pdf(paste(filename, "_error.pdf", sep = ""), height = 5, width = 7, pointsize = 12)
  plotErrorCurve(data, cols)
  dev.off() 
}

plotMappingBiasBenchmark <- function(data, cols, filename, min_count) {

  pdf(paste(filename, ".pdf", sep = ""), height = 4, width = 9, pointsize = 12)
  plotBiasCurve(data, cols, min_count)
  dev.off() 
}

plotMappingBiasBinomBenchmark <- function(data, cols, filename, alpha_thres, min_count) {
  
  pdf(paste(filename, ".pdf", sep = ""), height = 5, width = 9, pointsize = 12)
  plotBiasBinom(data, cols, alpha_thres, min_count)
  dev.off() 
}

plotMappingStatsBenchmark <- function(data, cols, filename) {

  pdf(paste(filename, ".pdf", sep = ""), height = 4, width = 4, pointsize = 12)
  plotStatsBarTwoLayer(data, cols)
  dev.off() 
}

plotMappingStatsBenchmarkWide <- function(data, cols, filename) {
  
  pdf(paste(filename, ".pdf", sep = ""), height = 4, width = 5, pointsize = 12)
  plotStatsBarTwoLayer(data, cols)
  dev.off() 
}

plotIsoSeqCorrelationBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    rename(Value = Corr) 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
  plotMapQCurve(data, cols, "Iso-Seq exon coverage correlation")
  dev.off() 
}

plotIsoSeqCoverageBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    rename(value_x = Coverage.est) %>%
    rename(value_y = Coverage.pb)
  
  pdf(paste(filename, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
  plotLogBins(data, cols, "Estimated exon coverage (log10 + 1)", "Iso-Seq exon coverage (log10 + 1)")
  dev.off() 
}

plotMappingComputeBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    rename(Value = Time)
  
  pdf(paste(filename, ".pdf", sep = ""), height = 4, width = 4, pointsize = 12)
  plotBar(data, cols, "Read pairs mapped per second")
  dev.off() 
}

plotMappingMemoryBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    rename(Value = Mem)
  
  pdf(paste(filename, ".pdf", sep = ""), height = 4, width = 4, pointsize = 12)
  plotBar(data, cols, "Maximum memory usage (GiB)")
  dev.off() 
}
