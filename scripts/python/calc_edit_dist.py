#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:34:51 2022

@author: Jordan
"""

import json
import sys

def edit_dist(aln):
    dist = 0
    for mapping in aln["path"]["mapping"]:
        for edit in mapping["edit"]:
            if "from_length" not in edit or edit["from_length"] == 0:
                if "to_length" in edit:
                    dist += int(edit["to_length"])
            elif "to_length" not in edit or edit["to_length"] == 0:
                if "from_length" in edit:
                    dist += int(edit["from_length"])
            elif "sequence" in edit and len(edit["sequence"]) != 0:
                dist += len(edit["sequence"])
    return dist

def cigar(aln):
    c = []
    for mapping in aln["path"]["mapping"]:
        for edit in mapping["edit"]:
            if "from_length" not in edit or edit["from_length"] == 0:
                if "to_length" in edit:
                    if len(c) != 0  and c[-1][1] == "I":
                        c[-1][0] += int(edit["to_length"])
                    else:
                        c.append([int(edit["to_length"]), "I"])
            elif "to_length" not in edit or edit["to_length"] == 0:
                if "from_length" in edit:
                    if len(c) != 0  and c[-1][1] == "D":
                        c[-1][0] += int(edit["from_length"])
                    else:
                        c.append([int(edit["from_length"]), "D"])
            elif "sequence" in edit and len(edit["sequence"]) != 0:
                if len(c) != 0  and c[-1][1] == "X":
                    c[-1][0] += int(edit["from_length"])
                else:
                    c.append([int(edit["from_length"]), "X"])
            else:
                if len(c) != 0  and c[-1][1] == "=":
                    c[-1][0] += int(edit["from_length"])
                else:
                    c.append([int(edit["from_length"]), "="])
                    
    return "".join("{}{}".format(l, op) for l, op in c)

if __name__ == "__main__":
    
    if len(sys.argv) == 1:
        json_in = sys.stdin
    else:
        assert(len(sys.argv) == 2)
        json_in = open(sys.argv[1])
        
    print("read_name\tcigar\tedit_dist\tscore")
    for line in json_in:
        aln = json.loads(line.strip())
        score = 0
        if "score" in aln:
            score = int(aln["score"])
        cig = cigar(aln)
        print("{}\t{}\t{}\t{}".format(aln["name"], cig, edit_dist(aln), score))