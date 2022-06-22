#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:58:38 2022

@author: Jordan
"""

import sys


if __name__ == "__main__":
    
    
    primary_mask = ~256
    
    for line in sys.stdin:
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith("@"):
            print(line.strip())
            continue
        
        tokens = line.strip().split("\t")
        
        tokens[1] = str(int(tokens[1]) & primary_mask)
        
        print("\t".join(tokens))