library(ggplot2)

setwd("/Users/Jordan/Documents/Research/Pangenomics/RNA/expr_analysis")
tab_fp = "rpvg_mpmap_r1_SRR12765650_1kg_all_af001_gencode100_unidi.top2.txt"
out_img = "rpvg_mpmap_r1_SRR12765650_1kg_all_af001_gencode100_unidi.png"
samp = "SRR12765650"
anc = "African American"
# tab_fp = "rpvg_mpmap_r1_SRR12765534_1kg_all_af001_gencode100_unidi.top2.txt"
# out_img = "rpvg_mpmap_r1_SRR12765534_1kg_all_af001_gencode100_unidi.png"
# samp = "SRR12765534"
# anc = "European American"

tab = read.table(tab_fp, header = 1)
min_reads = 1
tab = tab[tab$ReadCount >= min_reads,]

ord = order(tab$ReadCount)
rev_ord = rev(ord)
x = tab$ReadCount[ord]
prop_concentrated = function(thresh) rev(cumsum(tab$Top2ReadCount[rev_ord] >= thresh * tab$ReadCount[rev_ord]) / 1:nrow(tab))

#plot(x, prop_concentrated(.8), type = "l", ylim = c(0, 1), log = "x")

pre_plot_dat = data.frame(x)
colnames(plot_dat) = "x"
thresholds = c(.5, .6, .7, .8, .9, .99)
for (thresh in thresholds) {
    pre_plot_dat[[paste0("y", thresh)]] = prop_concentrated(thresh)
}

plot_dat = data.frame(rep(pre_plot_dat$x, length(thresholds)),
                      unlist(pre_plot_dat[2:ncol(pre_plot_dat)]),
                      rep(paste0(100 * thresholds, "%"), each = nrow(pre_plot_dat)))
colnames(plot_dat) = c("ReadCount", "Concentration", "Threshold")
plot_dat[["Threshold"]] = as.factor(plot_dat$Threshold)

plt_title = paste0(anc, " ancestry sample (", samp, ")")
colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
plt = (ggplot(plot_dat, aes(x = ReadCount, y = Concentration)) 
       + geom_line(aes(color = Threshold))
       + scale_x_continuous(trans='log10') 
       + ylim(0, 1)
       + scale_color_manual(values = colors)
       + ggtitle(plt_title)
       + ylab("Proportion of transcripts")
       + xlab("Minimum number of reads assigned to transcript") 
       + theme(plot.title = element_text(hjust = 0.5))
)
ggsave(out_img, plot = plt, dpi = 450, height = 5.5, width = 7.5, units = "in")
