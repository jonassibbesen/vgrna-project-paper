#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:19:09 2021

@author: Jordan
"""

import sys

if __name__ == "__main__":
    
    in_fp = sys.argv[1]
    
    curr_tx = None
    curr_clust = None
    total_ex = 0.0
    all_ex = None
    
    print("\t".join(["Name", "ClusterID", "ReadCount", "Top2ReadCount"]))
    
    def output_row():
        if curr_tx is not None:
            all_ex.sort()
            if len(all_ex) >= 2:
                top2 = all_ex[-1] + all_ex[-2]
            else:
                top2 = sum(all_ex)
            print("\t".join(map(str, [curr_tx, curr_clust, total_ex, top2])))
        
    with open(in_fp) as f:
        for line in f:
            if line.startswith("Name"):
                continue
            tokens = line.strip().split()
            tx = tokens[0].split("_")[0]
            clust = tokens[1]
            if tx != curr_tx:
                output_row()
                curr_tx = tx
                total_ex = 0.0
                all_ex = []
                curr_clust = clust
            
            ex = float(tokens[5])
            all_ex.append(ex)
            total_ex += ex
        output_row()