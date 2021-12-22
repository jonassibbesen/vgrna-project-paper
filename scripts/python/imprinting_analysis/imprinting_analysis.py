#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 12:12:25 2021

@author: Jordan
"""
import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
import tempfile
import gc
import re
import collections
import gzip
import bisect
import pickle
import itertools
import math


sns.set_style('whitegrid')

# make 2 maps:
# - from transcript ID to row numbers of corresponding haplotype specific transcripts
# - from cluster ID to corresponding transcript IDs
def row_dicts(tab):
    tx_rows = {}
    cluster_txs = {}
    for i in range(tab.shape[0]):
        tx_id = tab.Name.values[i].split("_")[0]
        clust_id = tab.ClusterID.values[i]
        if tx_id not in tx_rows:
            tx_rows[tx_id] = []
        tx_rows[tx_id].append(i)
        if clust_id not in cluster_txs:
            cluster_txs[clust_id] = set()
        cluster_txs[clust_id].add(tx_id)
    for clust_id in cluster_txs:
        cluster_txs[clust_id] = sorted(cluster_txs[clust_id])
    return tx_rows, cluster_txs

def gene_to_row_dict(tx_rows):
    gene_to_tx_rows = {}
    for tx_id in tx_id_to_gene:
        gene = tx_id_to_gene[tx_id]
        if gene not in gene_to_tx_rows:
            gene_to_tx_rows[gene] = []
        if tx_id in tx_rows:
            gene_to_tx_rows[gene].extend(tx_rows[tx_id])
    return gene_to_tx_rows

def parse_attr(attr):
    attrs = {}
    for t in attr.split(";"):
        tokens = t.strip().replace("\"", "").split()
        if len(tokens) == 0:
            continue
        tag, val = tokens
        attrs[tag] = val
    return attrs

def get_haplotypes(chrom, start, end, sample, genotypes):
    
    chrom_start = bisect.bisect_left(genotypes.CHROM.values, chrom)
    chrom_end = bisect.bisect_right(genotypes.CHROM.values, chrom)
    region_start = bisect.bisect_left(genotypes.POS.values, start, chrom_start, chrom_end)
    region_end = bisect.bisect_right(genotypes.POS.values, end, chrom_start, chrom_end)
    
    blocks = []
    for i in range(region_start, region_end):
        genotype = genotypes[sample].values[i]
        phased = "|" in genotype
        if len(blocks) == 0 or not phased:
            blocks.append({})
        
        al1, al2 = re.split("[\\|\\\\]", genotype)
        formatted_alleles = []
        for al in (al1, al2):
            fal = ""
            if al.isdigit():
                j = int(al)
                if j == 0:
                    fal = genotypes.REF.values[i]
                else:
                    fal = genotypes.ALT.values[i].split(",")[j - 1]
            formatted_alleles.append(fal)
        
        blocks[-1][genotypes.POS.values[i]] = tuple(formatted_alleles)
        
    return blocks


if __name__ == "__main__":
    
    assert(len(sys.argv) == 9)
        
    # gencode annotations
    gtf = sys.argv[1]
    # list of genes we're interested in
    focal_genes = sys.argv[2]
    # structured string in format SAMPLE1:rpvg_table1,SAMPLE2:rpvg_table2
    tab_string = sys.argv[3]
    # structured string in format SAMPLE1:sorted_gibbs_table1,SAMPLE2:sorted_gibbs_table2
    gibbs_string = sys.argv[4]
    # file constaining list of hst to variant files
    hst_variant_list = sys.argv[5]
    # file containing list of VCFs (probably reduced to these samples)
    vcf_list = sys.argv[6]
    # variants for the focal genes in one table
    variant_table = sys.argv[7] 
    # directory for output
    out_dir = sys.argv[8]
    
    
    tabs = []
    samples = []
    for tab_sample in tab_string.split(","):
        assert(":" in tab_sample)
        samp, tab = tab_sample.split(":")
        tabs.append(tab)
        samples.append(samp)
    gibbs_tabs = []
    gibbs_samples = []
    for tab_sample in gibbs_string.split(","):
        assert(":" in tab_sample)
        samp, tab = tab_sample.split(":")
        gibbs_tabs.append(tab)
        gibbs_samples.append(samp)
    assert(samples == gibbs_samples)
    
    assert(os.path.isdir(out_dir))
    assert(os.path.exists(gtf))
    assert(os.path.exists(focal_genes))
    assert(os.path.exists(vcf_list))
    for tab in tabs:
        assert(os.path.exists(tab))
    for tab in gibbs_tabs:
        assert(os.path.exists(tab))
    vcfs = []
    with open(vcf_list) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            vcf = line.strip()
            assert(os.path.exists(vcf))
            vcfs.append(vcf)
        
    # make a look table for the file name by chromosome
    hst_variant_files = {}
    with open(hst_variant_list) as f:
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8")
            fname = line.strip()
            with open(fname) as hst_f:
                #skip the header
                next(hst_f)
                hst_line = next(hst_f)
                if type(hst_line) == bytes:
                    hst_line = hst_line.decode("utf-8")
                hst_variant_files[hst_line.split()[0]] = fname
    
    
    
    tmpdir = tempfile.TemporaryDirectory()
    tmppref = tmpdir.name
    
    ###############
    
    focal_genes_set = set()
    for line in open(focal_genes):
        if type(line) == bytes:
            line = line.decode("utf-8") 
        focal_genes_set.add(line.strip().split()[0])
    
    ###############
    
    # load the GTF
    
    gencode = pd.read_csv(gtf, sep = "\t", header = None, skiprows = list(range(5)))
    gencode.columns = ["chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"]
    gencode['chr'] = gencode['chr'].apply(str)
    
    
    ###############
    
    print("loading gene annotations...", file = sys.stderr)
    
    # parse the GTF into useful indexes
    
    gene_coords = {}
    tx_models = {}
    tx_id_to_name = {}
    tx_id_to_gene = {}
    exonic_regions = {}
    
    for i in range(gencode.shape[0]):
        attrs = parse_attr(gencode.attr.values[i])
        gene = attrs["gene_id"]
        if gene not in tx_models:
            tx_models[gene] = {}
        chrom = gencode.chr.values[i]
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        if gene in tx_models:
            if gencode.type.values[i] == "gene":
                gene_coords[gene] = (chrom, gencode.start.values[i], gencode.end.values[i])
            elif gencode.type.values[i] == "exon":
                tx_id = attrs["transcript_id"]
                if tx_id not in tx_models[gene]:
                    tx_models[gene][tx_id] = []
                tx_models[gene][tx_id].append((chrom, gencode.start.values[i], gencode.end.values[i]))
                
                ###############
                
                tx_id_to_gene[tx_id] = gene
        
        ###############
        
        if "transcript_id" in attrs and "transcript_name" in attrs:
            tx_id_to_name[attrs["transcript_id"]] = attrs["transcript_name"]
        
        ###############
        
        if gencode.type.values[i] == "exon":
            if chrom not in exonic_regions:
                exonic_regions[chrom] = []
            exonic_regions[chrom].append([gencode.start.values[i], gencode.end.values[i]])
    
    ###############
    
    # reverse the transcript gene table
    
    gene_to_tx_ids = {}
    for tx_id in tx_id_to_gene:
        gene = tx_id_to_gene[tx_id]
        if gene not in gene_to_tx_ids:
            gene_to_tx_ids[gene] = []
        gene_to_tx_ids[gene].append(tx_id)
        
    ###############
    
    all_genes = sorted(gene_to_tx_ids)
    
    ###############
    
    # collapse the exonic regions that overlap
    
    for chrom in exonic_regions:
        i, j = 0, 0
        intervals = exonic_regions[chrom]
        intervals.sort()
        while j < len(intervals):
            if intervals[j][0] <= intervals[i][1]:
                intervals[i][1] = max(intervals[i][1], intervals[j][1])
            else:
                i += 1
                intervals[i] = intervals[j]
            j += 1
        while len(intervals) > i + 1:
            intervals.pop()
            
    ###############
    
    # this is a big table and we don't need it any more, clear it out
    
    del gencode
    gc.collect()
    
    ###############
    
    
    print("computing credible intervals...", file = sys.stderr)
    
    sample_tx_cred_intervals = {}
    for samp, tab in zip(gibbs_samples, gibbs_tabs):
        
        tx_cred_intervals = []
        sample_tx_cred_intervals[samp] = tx_cred_intervals
        
        def record_cred_interval(hst_exprs, credibility):
            if len(hst_exprs) == 0:
                return
            for hst1, hst2 in sorted(set(tuple(sorted(pair)) for pair in itertools.combinations(hst_exprs, 2))):
                ratios = []
                hst1_expr = hst_exprs[hst1]
                hst2_expr = hst_exprs[hst2]
                assert(len(hst1_expr) == len(hst2_expr))
                for i in range(len(hst1_expr)):
                    if hst1_expr[i] == 0.0 or hst2_expr[i] == 0.0:
                        # log ratio undefined if either is 0
                        continue
                    ratios.append(math.log(hst1_expr[i] / hst2_expr[i], 2.0))
                
                if len(ratios) == 0:
                    continue
                
                # find the credible interval
                ratios.sort()
                i1 = min(int(round(len(ratios) * (1.0 - credibility) / 2.0)), len(ratios) - 1)
                i2 = min(int(round(len(ratios) * (1.0 - (1.0 - credibility) / 2.0))), len(ratios) - 1)
                r1 = ratios[i1]
                r2 = ratios[i2]
                tx_cred_intervals.append((hst1, hst2, r1, r2))
                    
        # take either gzip or unzipped file
        f = None
        if tab.endswith(".gz"):
            f = gzip.open(tab)
        else:
            f = open(tab)
            
        # the credibility i'm using
        credibility = .9
    
        curr_tx = None
        hst_gibbs_exprs = None
        txs_seen = set()
        for line in f:
            if type(line) == bytes:
                line = line.decode("utf-8") 
            if line.startswith("Name"):
                # skip the header
                continue
            tokens = line.split()
            hst = tokens[0]
            tx = hst.split("_")[0]
            if tx != curr_tx:
                # were on to a new transcript, make sure we haven't seen it before
                assert(tx not in txs_seen)
                txs_seen.add(tx)
                
                if curr_tx is not None:
                    # record the ratios of the HSTs for the previous transcript
                    record_cred_interval(hst_gibbs_exprs, credibility)
                    
                # fresh data structures for this transcript
                curr_tx = tx
                hst_gibbs_exprs = {}
            
            # record the row of expression values
            hst_gibbs_exprs[hst] = [float(tokens[i]) for i in range(2, len(tokens))]
            
        if curr_tx is not None:
            # the final transcript
            record_cred_interval(hst_gibbs_exprs, credibility)
            
    
    sample_tx_cred_intervals_output = os.path.join(out_dir, "sample_tx_cred_intervals.pkl")
    with open(sample_tx_cred_intervals_output, "wb") as f:
        pickle.dump(sample_tx_cred_intervals, f)
    
    ###############
    
    
    print("loading genotypes...", file = sys.stderr)
    
    genotypes = pd.read_csv(variant_table, sep = "\t")
    genotypes['CHROM'] = genotypes['CHROM'].apply(str)
    genotypes.sort_values(["CHROM", "POS"], inplace = True)
    genotypes = genotypes.loc[np.invert(genotypes.duplicated()),:]
    
    
    
    #################
    
    print("loading HST variants...", file = sys.stderr)
    
    hst_variants = {}
    for hst_file in hst_variant_files.values():
        hst_table = pd.read_csv(hst_file, sep = "\t", header = 0)
        hst_table['Chrom'] = hst_table['Chrom'].apply(str)
        for i in range(hst_table.shape[0]):
            if type(hst_table.HSTs.values[i]) == float:
                # this seems to happen when the list of HSTs is empty
                continue
            hsts = hst_table.HSTs.values[i].split(",")
            for hst in hsts:
                tx = hst.split("_")[0]
                gene = tx_id_to_gene[tx]
                if not gene in focal_genes_set:
                    continue
                if not hst in hst_variants:
                    hst_variants[hst] = []
                var = (hst_table.Pos.values[i], hst_table.Allele.values[i])
                hst_variants[hst].append(var)
            

        del hst_table
        gc.collect()
    
    #################

    sample_higher_haplo_expr = {}
    sample_lower_haplo_expr = {}
    sample_informative_expr = {}
    sample_haplo_1_is_higher = {}
    sample_haplo_hsts = {}
    for i in range(len(tabs)):
        
        sample = samples[i]
        tab = tabs[i]
        
        print("computing haplotype expression for sample {}...".format(sample), file = sys.stderr)
        
        sample_expr = pd.read_csv(tab, sep = "\t")
        sample_tx_rows, sample_cluster_txs = row_dicts(sample_expr)
        
        higher_haplo_expr = {}
        lower_haplo_expr = {}
        informative_expr = {}
        haplo_1_is_higher = {}
        haplo_hsts = {}
        sample_higher_haplo_expr[sample] = higher_haplo_expr
        sample_lower_haplo_expr[sample] = lower_haplo_expr
        sample_informative_expr[sample] = informative_expr
        sample_haplo_1_is_higher[sample] = haplo_1_is_higher
        sample_haplo_hsts[sample] = haplo_hsts
        
        for gene in focal_genes_set:
            
            chrom, start, end = gene_coords[gene]
            blocks = get_haplotypes(chrom, start, end, sample, genotypes)
            
            if len(blocks) > 1:
                print("sample {} has {} phase blocks on gene {}, skipping".format(sample, len(blocks), gene), file = sys.stderr)
                continue
            
            block = blocks[0]
            
            if not gene in higher_haplo_expr:
                higher_haplo_expr[gene] = {}
                lower_haplo_expr[gene] = {}
                informative_expr[gene] = {}
            
            gene_higher_haplo_expr = higher_haplo_expr[gene]
            gene_lower_haplo_expr = lower_haplo_expr[gene]
            gene_informative_expr = informative_expr[gene]
            
            haplo_1_expr = {}
            haplo_2_expr = {}
            
            for tx_id in gene_to_tx_ids[gene]:
                haplo_1_expr[tx_id] = 0.0
                haplo_2_expr[tx_id] = 0.0
                total_informative_expr = 0.0
                haplo_hsts[tx_id] = [None, None]
                
                for i in sample_tx_rows[tx_id]:
                    ex = sample_expr.TPM.values[i]
                    hst = sample_expr.Name.values[i]
                    
                    match_1 = True
                    match_2 = True
                    for pos, allele in hst_variants[hst]:
                        hap_1, hap_2 = block[pos]
                        match_1 = match_1 and allele == hap_1
                        match_2 = match_2 and allele == hap_2
    
                        
                    if match_1 and not match_2:
                        haplo_hsts[tx_id][0] = hst
                        haplo_1_expr[tx_id] += ex
                    elif match_2 and not match_1:
                        haplo_hsts[tx_id][1] = hst
                        haplo_2_expr[tx_id] += ex
                    
                    if not (match_1 and match_2):
                        total_informative_expr += ex
                
                if not tx_id in gene_informative_expr:
                    gene_informative_expr[tx_id] = []
                gene_informative_expr[tx_id].append(total_informative_expr)
                        
                        
            if sum(haplo_1_expr.values()) > sum(haplo_2_expr.values()):
                higher = haplo_1_expr
                lower = haplo_2_expr
                haplo_1_is_higher[gene] = True
            else:
                lower = haplo_1_expr
                higher = haplo_2_expr
                haplo_1_is_higher[gene] = False
            
            for tx_id in higher:
                if not tx_id in gene_higher_haplo_expr:
                    gene_higher_haplo_expr[tx_id] = []
                    gene_lower_haplo_expr[tx_id] = []
                
                gene_higher_haplo_expr[tx_id].append(higher[tx_id])
                gene_lower_haplo_expr[tx_id].append(lower[tx_id])
    
    #################
    
    higher_haplo_output = os.path.join(out_dir, "sample_higher_haplo_expr.pkl")
    with open(higher_haplo_output, "wb") as f:
        pickle.dump(sample_higher_haplo_expr, f)
        
    lower_haplo_output = os.path.join(out_dir, "sample_lower_haplo_expr.pkl")
    with open(lower_haplo_output, "wb") as f:
        pickle.dump(sample_lower_haplo_expr, f)
        
    informative_output = os.path.join(out_dir, "sample_informative_expr.pkl")
    with open(informative_output, "wb") as f:
        pickle.dump(sample_informative_expr, f)
        
    which_haplo_output = os.path.join(out_dir, "sample_haplo_1_is_higher.pkl")
    with open(which_haplo_output, "wb") as f:
        pickle.dump(sample_haplo_1_is_higher, f)
        
    haplo_hsts_output = os.path.join(out_dir, "sample_haplo_hsts.pkl")
    with open(haplo_hsts_output, "wb") as f:
        pickle.dump(sample_haplo_hsts, f)
    
    ###############
    
    print("identifying heterozygous variants...", file = sys.stderr)
    
    inf = 2**62
    
    het_positions = {}
            
    for vcf in vcfs:
        with gzip.open(vcf) as f:
            samps = None
            for line in f:
                if type(line) == bytes:
                    line = line.decode("utf-8") 
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    samps = line.rstrip().split("\t")[9:]
                    for sample in samps:
                        if sample not in het_positions:
                            het_positions[sample] = set()
                else:
                    tokens = line.rstrip().split("\t")
                    assert(len(tokens) == len(samps) + 9)
                    chrom_exonic_regions = exonic_regions[tokens[0]]
                    chrom = tokens[0]
                    pos = int(tokens[1])
                    idx = bisect.bisect(chrom_exonic_regions, [pos, inf])
                    if idx == 0:
                        # before the first exon
                        continue
                    elif chrom_exonic_regions[idx - 1][1] < pos:
                        # in between exons
                        continue
                    for i in range(9, len(tokens)):
                        genotype = tokens[i]
                        samp = samps[i - 9]
                        if "|" in genotype or "\\" in genotype:
                            al1, al2 = re.split("[\\|\\\\]", genotype)
                            if al1 != al2:
                                het_positions[samp].add((chrom, pos))
        gc.collect()
    
    
    ###############
    
    all_gene_intervals = sorted((interval[0], interval[1], interval[2], gene) for gene, interval in gene_coords.items())
        
    sample_het_balance = {}
    for i in range(len(tabs)):
        
        tab = tabs[i]
        sample = samples[i]
        if sample not in sample_het_balance:
            sample_het_balance[sample] = {}
        het_balance = sample_het_balance[sample]
        
        print("computing balance for sample {}".format(sample), file = sys.stderr)
        
        buffer = collections.deque()
        prev_chrom = None
        tokens = None
        pos = None
        filesize = None
        hst_file = None
        
        gene_num = 0
        sample_expr = pd.read_csv(tab, sep = "\t")
        sample_tx_rows, sample_cluster_txs = row_dicts(sample_expr)
        
        for chrom, start, end, gene in all_gene_intervals:
    
            gene_num += 1
            if gene_num % 2500 == 0:
                print("processing gene {}".format(gene_num), file = sys.stderr)
                
            gene_hst_variants = {}
    
            if prev_chrom != chrom:
                # we've switched chromosomes to a new file
                if not chrom in hst_variant_files:
                    continue
                hst_table = hst_variant_files[chrom] 
                #print("starting chrom {}".format(chrom), file = sys.stderr)
                hst_file = open(hst_table)
                filesize = os.fstat(hst_file.fileno()).st_size
    
                # skip the header
                hst_file.readline()
                buffer.clear()
    
                tell = hst_file.tell()
                prev_pos = -1
                tokens = hst_file.readline().strip().split()
                var_chrom = tokens[0]
                pos = int(tokens[1])
                buffer.append((pos, tell))
    
    
    
            # advance through rows that are strictly before this gene
            while pos < start:
                tell = hst_file.tell()
                if tell == filesize:
                    break
                prev_pos = pos
                tokens = hst_file.readline().strip().split()
                pos = int(tokens[1])
                if pos != prev_pos:
                    buffer.append((pos, tell))
    
            # remove any part of the buffer before this gene
            while len(buffer) > 0:
                buf_pos = buffer[0][0]
                if buf_pos < start:
                    buffer.popleft()
                else:
                    break
    
            if len(buffer) > 0:
                # everything before the start has been removed, except the current row
                buf_pos, tell = buffer[0]
                if buf_pos < pos:
                    # this occurred strictly before the current row, so we need to seek
                    # backwards
    
                    # reset the part of the buffer to the right of where we're seeking to
                    while len(buffer) > 1:
                        buffer.pop()
                    hst_file.seek(tell)
    
                    tokens = hst_file.readline().strip().split()
                    pos = int(tokens[1])
    
                    
            hst_vars = {}
            # iterate over rows in the gene
            while pos <= end:
                
                if len(tokens) >= 5:
                    allele = tokens[3]
                    pos = int(tokens[1])
                    hsts = tokens[4].split(",")
                    for hst in hsts:
                        if hst not in hst_vars:
                            hst_vars[hst] = []
                        hst_vars[hst].append((pos, allele))
    
                tell = hst_file.tell()
                if tell == filesize:
                    # we hit the end of the file
                    break
                prev_pos = pos
                tokens = hst_file.readline().strip().split()
                pos = int(tokens[1])
                if pos != prev_pos:
                    # this is the first row we've seen with this position, remember
                    # it in the buffer
                    buffer.append((pos, tell))
    
    
            prev_chrom = chrom
            
            if gene not in het_balance:
                het_balance[gene] = []
            
            var_expr = {}
            if gene not in gene_to_tx_ids:
                continue
            for tx_id in gene_to_tx_ids[gene]:
                #print("looking at expression for tx " + tx_id, file = sys.stderr)
                if tx_id not in sample_tx_rows:
                    continue
                for i in sample_tx_rows[tx_id]:
                    ex = sample_expr.TPM.values[i]
                    if ex == 0.0:
                        continue
                    hst = sample_expr.Name.values[i]
                    #print("\thst " + hst + " has positive expression " + str(ex), file = sys.stderr)
                    if hst not in hst_vars:
                        # must not overlap any variants
                        continue
                    for var in hst_vars[hst]:
                        if var not in var_expr:
                            var_expr[var] = 0.0
                        var_expr[var] += ex
        
            alleles = {}
            for pos, allele in var_expr:
                if pos not in alleles:
                    alleles[pos] = []
                alleles[pos].append(allele)
    
                
            for pos in alleles:
                if (chrom, pos) not in het_positions[sample]:
                    continue
                #print("looking at expression for pos " + chrom + " " + str(pos), file = sys.stderr)
                total_expr = sum(var_expr[(pos, allele)] for allele in alleles[pos])
                highest_expr = max(var_expr[(pos, allele)] for allele in alleles[pos])
                #print("highest expr " + str(highest_expr) + ", total " + str(total_expr), file = sys.stderr)
                
                het_balance[gene].append((highest_expr, total_expr))
                
        del sample_expr
        del sample_tx_rows
        del sample_cluster_txs
        gc.collect()
    
    
    #################
    
    balance_output = os.path.join(out_dir, "sample_het_balance.pkl")
    with open(balance_output, "wb") as f:
        pickle.dump(sample_het_balance, f)
        
    tx_models_output = os.path.join(out_dir, "tx_models.pkl")
    with open(tx_models_output, "wb") as f:
        pickle.dump(tx_models, f)
        
    tx_id_to_name_output = os.path.join(out_dir, "tx_id_to_name.pkl")
    with open(tx_id_to_name_output, "wb") as f:
        pickle.dump(tx_id_to_name, f)
    
    
