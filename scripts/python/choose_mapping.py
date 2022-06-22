#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:19:54 2022

@author: Jordan
"""

import sys
import random

SECONDARY = 256
SUPPLEMENTARY = 2048

def to_str(s):
    if type(s) == bytes:
        return s.decode("utf-8")
    else:
        return s

# print the merged headers and 
def process_headers(f1, f2):
    
    # lines above the SQ lines
    line1 = to_str(next(f1))
    while not line1.startswith("@SQ"):
        print(line1, end = "")
        line1 = to_str(next(f1))
    line2 = to_str(next(f2))
    while not line2.startswith("@SQ"):
        line2 = to_str(next(f2))
    
    # SQ lines
    while line1.startswith("@SQ"):
        print(line1, end = "")
        line1 = to_str(next(f1))
    while line2.startswith("@SQ"):
        if not "KI" in line2 and not "GL" in line2:
            print(line2, end = "")
        line2 = to_str(next(f2))
        
    # lines after SQ lines
    while line1 is not None and line1.startswith("@"):
        print(line1, end = "")
        line1 = to_str(next(f1))
    
    while line2 is not None and line2.startswith("@"):
        line2 = to_str(next(f2))
        
    return line1.strip().split("\t"), line2.strip().split("\t")
    

# get next sam record, skipping secondary and supplementary alignments
def get_next(f):
    try:
        r = next(f)
    except StopIteration:
        return None
    else:
        r = to_str(r)
        tokens = r.strip().split("\t")
        flags = int(tokens[1])
        while flags & SECONDARY or flags & SUPPLEMENTARY:
            try:
                r = next(f)
            except StopIteration:
                return None
            r = to_str(r)
            tokens = r.strip().split("\t")
            flags = int(tokens[1])
            
    return tokens
            
def score(tokens):
    s = 0
    if len(tokens) >= 12:
        for tag in tokens[11:]:
            if tag.startswith("AS"):
                s = int(tag.split(":")[2])
                break
    return s

if __name__ == "__main__":
    
    sam1 = sys.argv[1]
    sam2 = sys.argv[2]
    
    f1 = open(sam1)
    f2 = open(sam2)
    
    random.seed(982316498213)
    
    tokens1, tokens2 = process_headers(f1, f2)
    
    while tokens1 is not None and tokens2 is not None:
        
        # qnames should match
        assert(tokens1[0] == tokens2[0])
        
        score1 = score(tokens1)
        score2 = score(tokens2)
        
        emit = None
        if score1 > score2:
            emit = tokens1
        elif score1 < score2:
            emit = tokens2
        else:
            emit = random.choice([tokens1, tokens2])
            
        print("\t".join(emit))
        
        tokens1 = get_next(f1)
        tokens2 = get_next(f2)
    
    