#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:42:38 2021

@author: Jordan
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import math
import random
import pickle
import bisect
import os

sns.set_style('whitegrid')

def bin_values(xs, bins):
    assert(list(bins) == sorted(bins))
    assert(len(bins) > 1)
    counts = [0 for b in range(len(bins) - 1)]
    for x in xs:
        i = bisect.bisect(bins, x)
        if i == 0:
            i += 1
        if i == len(bins):
            i -= 1
        counts[i - 1] += 1
    return counts

if __name__ == "__main__":
    
    pickle_dir = sys.argv[1]
    focal_genes = sys.argv[2]
    plot_dir = sys.argv[3]
    do_genes = True
    if len(sys.argv) == 5:
        do_genes = bool(int(sys.argv[4]))
    
    imprinted_genes = set()
    with open(focal_genes) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8") 
            imprinted_genes.add(line.strip().split()[0])
    
    with open(os.path.join(pickle_dir, "sample_het_balance.pkl"), "rb") as f:
        sample_het_balance = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_higher_haplo_expr.pkl"), "rb") as f:
        sample_higher_haplo_expr = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_lower_haplo_expr.pkl"), "rb") as f:
        sample_lower_haplo_expr = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_informative_expr.pkl"), "rb") as f:
        sample_informative_expr = pickle.load(f)
    with open(os.path.join(pickle_dir, "tx_models.pkl"), "rb") as f:
        tx_models = pickle.load(f)
    with open(os.path.join(pickle_dir, "tx_id_to_name.pkl"), "rb") as f:
        tx_id_to_name = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_haplo_hsts.pkl"), "rb") as f:
        sample_haplo_hsts = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_haplo_1_is_higher.pkl"), "rb") as f:
        sample_haplo_1_is_higher = pickle.load(f)
    with open(os.path.join(pickle_dir, "sample_tx_cred_intervals.pkl"), "rb") as f:
        sample_tx_cred_intervals = pickle.load(f)
        
        
        
        
    samples = sorted(sample_het_balance)
    
    min_tpm = 1.0
    sample_balance_all = {}
    sample_balance_imprinted = {}
    for sample in sample_het_balance:
        het_balance = sample_het_balance[sample]
        sample_balance_all[sample] = []
        sample_balance_imprinted[sample] = []
        for gene in het_balance:
            for max_ex, tot_ex in het_balance[gene]:
                if tot_ex < 1.0:
                    continue
                balance = max_ex / tot_ex
                sample_balance_all[sample].append(balance)
                if gene in imprinted_genes:
                    sample_balance_imprinted[sample].append(balance)
        
    resolution = 300
    
    for sample in samples:
        imprinted_label = "Top 20 imprinted genes\nfrom Zink, et al. (2018)"
        all_label = "All genes"
        imprinted_color = "lightblue"
        all_color = "lightcoral"
        imprinted_label_color = "steelblue"
        all_label_color = "firebrick"
        border_color = "darkgray"
        bar_linewidth = .00
        
        balance_imprinted = sample_balance_imprinted[sample]
        balance_all = sample_balance_all[sample]
        
        bins = np.linspace(.4, 1, 25)
        balance_bins_imprinted = [bin_cnt / len(balance_imprinted) for bin_cnt in bin_values(balance_imprinted, bins)]
        balance_bins_all = [-bin_cnt / len(balance_all) for bin_cnt in bin_values(balance_all, bins)]
        
        x_max = 1.05 * max(max(balance_bins_imprinted), -min(balance_bins_all))
        x_min = -x_max
        y_max = 1.01
        y_min = 0.465
        
        bin_centers = [(bins[i] + bins[i + 1]) / 2.0 for i in range(len(bins) - 1)]
        
        delta = (bins[-1] - bins[0]) / (len(bins) - 1) - bar_linewidth / 2.0
        fig, ax = plt.subplots(1, 1, dpi = resolution)
        ax.barh(bin_centers, balance_bins_all, delta, align = "center", color = all_color, 
                label = all_label, linewidth = bar_linewidth)
        ax.barh(bin_centers, balance_bins_imprinted, delta, align = "center", color = imprinted_color, 
                label = imprinted_label, linewidth = bar_linewidth)
        ax.plot([0, 0], [y_max, y_min], border_color, linewidth=1)
        ax.tick_params(axis='x', colors=border_color) 
        
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.grid(False, axis = "x")
        
        all_tick_step = 4000
        imprinted_tick_step = 5
        
        x_tick = []
        x_tick_label = []
        x_tick_label_col = []
        for x in range(int(math.ceil(x_min * len(balance_all) / all_tick_step)) * all_tick_step, 0, all_tick_step):
            x_tick.append(x / len(balance_all))
            x_tick_label.append(str(-x))
            x_tick_label_col.append(all_label_color)
        x_tick.append(0.0)
        x_tick_label.append("0")
        x_tick_label_col.append("black")
        for x in range(imprinted_tick_step, int(math.ceil(x_max * len(balance_imprinted))), imprinted_tick_step):
            x_tick.append(x / len(balance_imprinted))
            x_tick_label.append(str(x))
            x_tick_label_col.append(imprinted_label_color)
        
        ax.xaxis.set_label_position('top') 
        ax.xaxis.tick_top()
        ax.set_xticks(x_tick, minor=False)
        ax.set_xticklabels(x_tick_label, minor=False)
        for lab, color in zip(ax.get_xticklabels(), x_tick_label_col):
            lab.set_color(color)
        
        for spine in list(ax.spines.values()):
            spine.set_edgecolor(border_color)
            
        ax.spines['right'].set_visible(False)
        #ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        ax.set(xlabel='Number of heterozygous variants', ylabel='Proportion of expression')
        
        #plt.legend(bbox_to_anchor=(-0.02,1.015), loc="upper left", frameon=False)
        #plt.legend(bbox_to_anchor=(1.015,.035), loc="lower right", frameon=False)
        plt.legend(bbox_to_anchor=(0,.81), loc="upper left", frameon=False)
        #plt.legend(bbox_to_anchor=(.57,.822), loc="upper left", frameon=False)
        plt.savefig(os.path.join(plot_dir, sample + "_bidirectional_variant_expression.png"), 
                    bbox_inches="tight")
        
    #  decrease default font size of 10
    plt.rcParams.update({'font.size': 8.5})
        
    if do_genes:
        for sample in samples:
            
            higher_haplo_expr = sample_higher_haplo_expr[sample]
            lower_haplo_expr = sample_lower_haplo_expr[sample]
            informative_expr = sample_informative_expr[sample]
            haplo_hsts = sample_haplo_hsts[sample]
            haplo_1_is_higher = sample_haplo_1_is_higher[sample]
            tx_cred_intervals = sample_tx_cred_intervals[sample]
            
            hst_intervals = {}
            for i in range(len(tx_cred_intervals)):
                hst1, hst2, r1, r2 = tx_cred_intervals[i]
                if not hst1 in hst_intervals:
                    hst_intervals[hst1] = []
                if not hst2 in hst_intervals:
                    hst_intervals[hst2] = []
                hst_intervals[hst1].append(i)
                hst_intervals[hst2].append(i)
            
            
            for gene in imprinted_genes:
                
                first_tx_name = tx_id_to_name[next(iter(tx_models[gene]))]
                gene_common_name = first_tx_name[:first_tx_name.rfind("-")]
                exon_thickness = .3
                intron_thickness = .075
                bars_fraction = 2.0 / 3.0
                expr_bars_fraction = .4
                shoulder = .25
                absolute = False
                min_tpm = .1
                scatter_cols = ["tab:red", "tab:purple", "tab:olive"]
                bar_cols = ["#882255", "#44AA99", "#DDCC77"]
                expr_bar_col = "gray"
                
                hap_1_higher = haplo_1_is_higher[gene]
                tx_higher_haplo_expr = higher_haplo_expr[gene]
                tx_lower_haplo_expr = lower_haplo_expr[gene]
                tx_total_informative_expr = informative_expr[gene]
                
                tx_max_expr = max(sum(v) for v in tx_total_informative_expr.values())
                
                txs = sorted([tx for tx in tx_total_informative_expr if sum(tx_total_informative_expr[tx]) >= min_tpm], 
                             key = lambda tx: sum(tx_total_informative_expr[tx]))
                
                tx_names = [tx_id_to_name[tx_id] for tx_id in txs]
                
#                f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize = (8, 4), dpi = resolution,
#                                                 gridspec_kw={'width_ratios': [1, 1, 1]})
                f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize = (11, 4), dpi = resolution,
                                                     gridspec_kw={'width_ratios': [1, 1, 1, 1]})
                
                chrom = None
                
                for j in range(len(txs)):
                    tx = txs[j]
                    tx_model = tx_models[gene][tx]
                    alpha = 1.0
                #     if sum(tx_total_informative_expr[tx]) != 0.0:
                #         alpha = sum(tx_total_informative_expr[tx]) / tx_max_expr
                #     else:
                #         alpha = 0.0
                    chrom = tx_model[0][0]
                    for i in range(len(tx_model)):
                        exon = patches.Rectangle((tx_model[i][1], j - .5 * exon_thickness), 
                                                 tx_model[i][2] - tx_model[i][1], exon_thickness, 
                                                 linewidth=0, edgecolor=None, facecolor='k', alpha = alpha)
                
                        ax1.add_patch(exon)
                        if i != 0:
                            if tx_model[i - 1][2] < tx_model[i][1]:
                                intron_begin = tx_model[i - 1][2]
                                intron_len = tx_model[i][1] - tx_model[i - 1][2]
                            else:
                                intron_begin = tx_model[i][2]
                                intron_len = tx_model[i - 1][1] - tx_model[i][2]
                            intron = patches.Rectangle((intron_begin, j - .5 * intron_thickness), 
                                                       intron_len, intron_thickness, 
                                                       linewidth=0, edgecolor=None, facecolor='k', alpha = alpha)
                            ax1.add_patch(intron) 
                
                xmin, xmax = 2**60, -1
                for tx in tx_models[gene]:
                    if tx not in txs:
                        continue
                    for chrom, start, end in tx_models[gene][tx]:
                        xmin = min(xmin, start)
                        xmax = max(xmax, end)  
                
                center = (xmax + xmin) / 2.0
                radius = (xmax - xmin) / 2.0
                xmax = center + (1.0 + shoulder) * radius
                xmin = center - (1.0 + shoulder) * radius
                
                ax1.set_xlim(xmin, xmax)
                
                y_pos = []
                for i in range(3):
                    y_pos.append([j - bars_fraction / 2.0 + bars_fraction / 6.0 + i * bars_fraction / 3.0 for j in range(len(txs))])
                
                x_highest = []
                x_second_highest = []
                x_rest = []
                x_scatter = []
                y_scatter = []
                c_scatter = []
                tx_tpm = []
                for i in range(len(txs)):
                    tx = txs[i]
                    tx_tpm.append(sum(tx_total_informative_expr[tx]))
                    if sum(tx_total_informative_expr[tx]) == 0.0 or tx not in tx_higher_haplo_expr:
                        x_highest.append(0.0)
                    else:
                        if absolute:
                            x_highest.append(sum(tx_higher_haplo_expr[tx]) / len(tx_higher_haplo_expr[tx]))
                        else:
                            x_highest.append(sum(tx_higher_haplo_expr[tx]) / sum(tx_total_informative_expr[tx]))
                        for j in range(len(tx_higher_haplo_expr[tx])):
                            haplo_ex = tx_higher_haplo_expr[tx][j]
                            inform_ex = tx_total_informative_expr[tx][j]
                            if absolute:
                                x_scatter.append(haplo_ex)
                            else:
                                if inform_ex != 0.0:
                                    x_scatter.append(haplo_ex / inform_ex)
                                else:
                                    x_scatter.append(0.0)
                            y_scatter.append(y_pos[2][i] + random.uniform(-bars_fraction / 12.0, bars_fraction / 12.0))
                            c_scatter.append(scatter_cols[j])
                for i in range(len(txs)):
                    tx = txs[i]
                    if sum(tx_total_informative_expr[tx]) == 0.0 or tx not in tx_lower_haplo_expr:
                        x_second_highest.append(0.0)
                    else:
                        if absolute:
                            x_second_highest.append(sum(tx_lower_haplo_expr[tx]) / len(tx_lower_haplo_expr[tx]))
                        else:
                            x_second_highest.append(sum(tx_lower_haplo_expr[tx]) / sum(tx_total_informative_expr[tx]))
                        for j in range(len(tx_higher_haplo_expr[tx])):
                            haplo_ex = tx_lower_haplo_expr[tx][j]
                            inform_ex = tx_total_informative_expr[tx][j]
                            if absolute:
                                x_scatter.append(haplo_ex)
                            else:
                                if inform_ex != 0.0:
                                    x_scatter.append(haplo_ex / inform_ex)
                                else:
                                    x_scatter.append(0.0)
                            y_scatter.append(y_pos[1][i] + random.uniform(-bars_fraction / 12.0, bars_fraction / 12.0))
                            c_scatter.append(scatter_cols[j])
                for i in range(len(txs)):
                    tx = txs[i]
                    if sum(tx_total_informative_expr[tx]) == 0.0:
                        x_rest.append(0.0)
                    else:
                        if absolute:
                            x_rest.append((sum(tx_total_informative_expr[tx]) - sum(tx_lower_haplo_expr[tx] + tx_higher_haplo_expr[tx])) / len(tx_total_informative_expr[tx]))
                        else:
                            x_rest.append(1.0 - sum(tx_lower_haplo_expr[tx] + tx_higher_haplo_expr[tx]) / sum(tx_total_informative_expr[tx]))
                        for j in range(len(tx_higher_haplo_expr[tx])):
                            haplo_ex_1 = tx_higher_haplo_expr[tx][j]
                            haplo_ex_2 = tx_lower_haplo_expr[tx][j]
                            inform_ex = tx_total_informative_expr[tx][j]
                            if absolute:
                                x_scatter.append(inform_ex - haplo_ex_1 - haplo_ex_2)
                            else:
                                if inform_ex != 0.0:
                                    x_scatter.append(1.0 - (haplo_ex_1 + haplo_ex_2) / inform_ex)
                                else:
                                    x_scatter.append(0.0)
                            y_scatter.append(y_pos[0][i] + random.uniform(-bars_fraction / 12.0, bars_fraction / 12.0))
                            c_scatter.append(scatter_cols[j])
                            
                # darker line at the zero mark
                ax4.plot([0, 0], [-1, len(txs) + 1], c = "darkgray", linewidth = .75)
                            
                interval_max = 0
                for i in range(len(txs)):
                    tx = txs[i]
                    
                    hst1, hst2 = haplo_hsts[tx]
                    if hst1 is None or hst2 is None:
                        print("credibility interval for transcript {} of gene {} skipped for no haplo HSTs".format(tx, gene), file = sys.stderr)
                        continue
                    if hap_1_higher:
                        hst_higher, hst_lower = hst1, hst2
                    else:
                        hst_higher, hst_lower = hst2, hst1
                        
                    r1, r2 = None, None
                    if hst_higher not in hst_intervals:
                        print("credibility interval for transcript {} of gene {} skipped for no sampled iterations".format(tx, gene), file = sys.stderr)
                        continue
                    for j in hst_intervals[hst_higher]:
                        if tx_cred_intervals[j][1] == hst_lower:
                            r1, r2 = tx_cred_intervals[j][2], tx_cred_intervals[j][3]
                            break
                        elif tx_cred_intervals[j][0] == hst_lower:
                            r1, r2 = -tx_cred_intervals[j][3], -tx_cred_intervals[j][2]
                            break
                        
                    if r1 == None:
                        print("credibility interval for transcript {} of gene {} skipped".format(tx, gene), file = sys.stderr)
                        continue
                    
                    interval_col = "black"
                    ax4.plot([r1, r2], [y_pos[1][i], y_pos[1][i]], c = interval_col)
                    ax4.plot([r1, r1], [(y_pos[1][i] + y_pos[2][i]) / 2.0, (y_pos[1][i] + y_pos[0][i]) / 2.0], c = interval_col)
                    ax4.plot([r2, r2], [(y_pos[1][i] + y_pos[2][i]) / 2.0, (y_pos[1][i] + y_pos[0][i]) / 2.0], c = interval_col)
                    
                    
                    interval_max = max(interval_max, r2, -r1)

                interval_max *= 1.05
                ax4.set_xlim(-interval_max, interval_max)
                
                ratio_tick_max = int(interval_max / math.log(10, 2.0))
                    
                ratio_ticks = []
                ratio_tick_labels = []
                for i in range(-ratio_tick_max, ratio_tick_max + 1):
                    ratio_ticks.append(i * math.log(10, 2.0))
                    if i <= 0:
                        ratio_tick_labels.append("1:" + str(10**(-i)))
                    else:
                        ratio_tick_labels.append(str(10**(i)) + ":1")
                ax4.xaxis.set_ticks(ratio_ticks)
                ax4.xaxis.set_ticklabels(ratio_tick_labels)
                    
                plt.setp(ax1.xaxis.get_majorticklabels(), rotation=0)
                
                #ax1.yaxis.set_visible(False)
                ax1.grid(which='major', axis='y')
                ax1.yaxis.tick_left()
                ax1.yaxis.set_ticks(y_pos[1])
                ax1.yaxis.set_ticklabels(tx_names)
                ax1.tick_params(axis='y', length=0.0)
                ax1.ticklabel_format(axis = "x", style = "plain")
                    
                ax2.yaxis.set_visible(False)
                ax3.yaxis.set_visible(False)
                ax4.yaxis.set_visible(False)
                
                ax1.set_ylim(-exon_thickness, len(txs) - 1 + exon_thickness)
                
                ax2.barh(y_pos[1], tx_tpm, expr_bars_fraction, align = "center", color = expr_bar_col,
                         linewidth = 0.0)
                
                bar_fraction_shrinkage = .95
                
                ax3.barh(y_pos[2], x_highest, bar_fraction_shrinkage * bars_fraction / 3.0, color = bar_cols[0], align='center', label = "Haplotype 1",
                         linewidth = 0.0)    
                ax3.barh(y_pos[1], x_second_highest, bar_fraction_shrinkage * bars_fraction / 3.0, color = bar_cols[1], align='center', label = "Haplotype 2",
                         linewidth = 0.0)    
                ax3.barh(y_pos[0], x_rest, bar_fraction_shrinkage * bars_fraction / 3.0, color = bar_cols[2], align='center', label='Incorrect haplotypes',
                         linewidth = 0.0)
                #     ax2.scatter(x_scatter, y_scatter, alpha = .75, c = c_scatter, s = .5, zorder = 100)
                if not absolute:
                    ax3.set_xlim(-.01, 1.01)
                else:
                    ax3.set_xlim(-.01 * max(x_scatter), 1.01 * max(x_scatter))
                    
                ax1.spines['top'].set_visible(False)
                ax1.spines['right'].set_visible(False)
                ax1.spines['bottom'].set_visible(False)
                ax1.spines['left'].set_visible(False)
                ax2.spines['top'].set_visible(False)
                ax2.spines['right'].set_visible(False)
                ax2.spines['bottom'].set_visible(False)
                ax3.spines['top'].set_visible(False)
                ax3.spines['right'].set_visible(False)
                ax3.spines['bottom'].set_visible(False)
                ax3.spines['left'].set_visible(False)
                ax4.spines['top'].set_visible(False)
                ax4.spines['right'].set_visible(False)
                ax4.spines['bottom'].set_visible(False)
                ax4.spines['left'].set_visible(False)
                
                ax1.set(xlabel='Chromosome {} position'.format(chrom), ylabel='')
                ax1.xaxis.set_label_coords(.333, -0.1)
                ax2.set(xlabel='Total expression (TPM)', ylabel='')
                ax3.set(xlabel='Proportion of expression', ylabel='')
                ax4.set(xlabel='Haplotype ratio\n(90% credible interval)', ylabel='')
                
                bar_labels = ["Haplotype 1", "Haplotype 2", "Incorrect haplotypes"]
                leg1 = plt.legend(handles = [patches.Patch(color = col, label = lab) for col, lab in zip(bar_cols, bar_labels)], 
                                  bbox_to_anchor=(-1.375,-.2), loc="center", frameon=False, ncol = len(bar_cols))
                ax4.add_artist(leg1)
                # leg2 = plt.legend(handles = [patches.Patch(color='black', alpha = k / 4.0, label="{:.2f}".format(tx_max_expr * k / 4.0)) for k in range(1, 5)],
                #                bbox_to_anchor=(1.04,0), loc="lower left", title="Avg. expression (TPM)")
                # ax3.add_artist(leg2)
                #f.suptitle("Haplotype expression for isoforms of " + gene_common_name)
                f.savefig(os.path.join(plot_dir, sample + "_" + gene_common_name + "_isoform_expression.png"), bbox_inches="tight")
    
    