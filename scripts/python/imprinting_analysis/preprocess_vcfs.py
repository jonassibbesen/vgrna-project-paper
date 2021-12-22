#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 13:07:48 2021

@author: Jordan
"""

import sys
import subprocess
import os
import shutil
import tempfile
import gc
import pandas as pd

bcftools = ""
tabix = ""
bgzip = ""

def parse_attr(attr):
    attrs = {}
    for t in attr.split(";"):
        tokens = t.strip().replace("\"", "").split()
        if len(tokens) == 0:
            continue
        tag, val = tokens
        attrs[tag] = val
    return attrs

def print_help():
    print("usage: ./preprocess_vcfs.py <GTF gene annotation> <table of focal genes> <list of VCF files> <list of samples to include> <dir for VCF output>", file = sys.stderr)
    print("prints VCF without header lines for selected samples and genes", file = sys.stderr)
    sys.exit()

if __name__ == "__main__":
    
    if len(sys.argv) != 6:
        print_help()
    
    gtf =  sys.argv[1]
    focal_genes =  sys.argv[2]
    vcf_list = sys.argv[3]
    sample_list = sys.argv[4]
    vcf_output_dir= sys.argv[5]
    
    if not os.path.isfile(vcf_list) or not os.path.isfile(sample_list) or not os.path.isfile(gtf) or not os.path.isfile(focal_genes) or not os.path.isdir(vcf_output_dir):
        print_help()

    if bcftools == "":
        bcftools = shutil.which("bcftools")
    if bgzip == "":
        bgzip = shutil.which("bgzip")
    if tabix == "":
        tabix = shutil.which("tabix")

    if bcftools is None or bgzip is None or tabix is None:
        print("must have bcftools, bgzip, and tabix installed", file = sys.stderr)
        sys.exit()
    
    tmpdir = tempfile.TemporaryDirectory()
    tmppref = tmpdir.name
    
    samples = []
    vcfs = []
    with open(vcf_list) as f:
        for line in f:
            vcfs.append(line.rstrip())
    with open(sample_list) as f:
        for line in f:
            samples.append(line.rstrip())
            
    ###############
    
    genes = set()
    for line in open(focal_genes):
        if type(line) == bytes:
            line = line.decode("utf-8") 
        genes.add(line.strip().split()[0])
    
    ###############
        
    print("extracting samples from VCFs...", file = sys.stderr)
    
    # extract the samples we want into new vcfs
    procs = []
    sample_vcf_names = []
    for i in range(len(vcfs)):
        vcf = vcfs[i]
        sample_vcf_name = os.path.join(vcf_output_dir, "sampled_vcf_" + str(i) + ".vcf")
        sample_vcf_names.append(sample_vcf_name)
        proc = subprocess.Popen([bcftools, "view", "-s", ",".join(samples), vcf], stdout = open(sample_vcf_name, "wb"))
        procs.append(proc)
        
    # barrier sync
    for proc in procs:
        proc.communicate()
    procs.clear()
    
    print("indexing VCFs for region extraction...", file = sys.stderr)
    
    # bgzip the sample vcfs
    for vcf in sample_vcf_names:
        proc = subprocess.Popen([bgzip, vcf])
        procs.append(proc)
        
    # barrier sync
    for proc in procs:
        proc.communicate()
    procs.clear()
        
    # bgzip the sample vcfs
    for vcf in sample_vcf_names:
        zipped_vcf = vcf + ".gz"
        proc = subprocess.Popen([tabix, "-p", "vcf", zipped_vcf])
        procs.append(proc)
        
    # barrier sync
    for proc in procs:
        proc.communicate()
    procs.clear()
    
    ###############
    
    # load the GTF
    
    print("loading annotations...", file = sys.stderr)
    
    gencode = pd.read_csv(gtf, sep = "\t", header = None, skiprows = list(range(5)))
    gencode.columns = ["chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"]
    gencode['chr'] = gencode['chr'].apply(str)
    
    ###############
    
    print("identifying gene regions...", file = sys.stderr)
    
    gene_coords = {}
    
    for i in range(gencode.shape[0]):
        attrs = parse_attr(gencode.attr.values[i])
        gene = attrs["gene_id"]
        if gene in genes:
            if gencode.type.values[i] == "gene":
                chrom = gencode.chr.values[i]
                if chrom.startswith("chr"):
                    chrom = chrom[3:]
                gene_coords[gene] = (chrom, gencode.start.values[i], gencode.end.values[i])
           
    ###############
    
    # this is a big table and we don't need it any more, clear it out
    
    del gencode
    gc.collect()
    
    ###############
    
    gene_bed_name = os.path.join(tmppref, "focal_genes.bed")
    with open(gene_bed_name, "w") as f:
        for gene in sorted(gene_coords, key = lambda gene: gene_coords[gene]):
            print("\t".join(str(v) for v in gene_coords[gene]), file = f)
    
    print("extracting gene regions...", file = sys.stderr)
    
    # subset by the vcfs by the focal genes
    region_vcfs = []
    for i in range(len(sample_vcf_names)):
        zipped_vcf = sample_vcf_names[i] + ".gz"
        region_vcf = os.path.join(tmppref, "sample_region_" + str(i) + ".vcf")
        region_vcfs.append(region_vcf)
        proc = subprocess.Popen([tabix, "-R", gene_bed_name, zipped_vcf], 
                                stdout = open(region_vcf, "wb"))
        procs.append(proc)
        
    # barrier sync
    for proc in procs:
        proc.communicate()
    procs.clear()
    
    print("writing genotype table to STDOUT...", file = sys.stderr)
        
    # get the header (ugly hack but it works)
    header_file = os.path.join(tmppref, "header.txt")
    proc = subprocess.Popen([bcftools, "view", "-h", sample_vcf_names[0] + ".gz"], 
                            stdout = open(header_file, "wb"))
    proc.communicate()
    with open(header_file, "rb") as f:
        for line in f:
            if not line.startswith(b"##"):
                sys.stdout.buffer.write(line[1:])
    
    # write out a table to stdout
    for i in range(len(region_vcfs)):
        with open(region_vcfs[i], "rb") as f:
            for line in f:
                if line.startswith(b"#"):
                    continue
                sys.stdout.buffer.write(line)
                    
