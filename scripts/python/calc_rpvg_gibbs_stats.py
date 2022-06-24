
'''
calc_rpvg_gibbs_stats.py 
Calculates rpvg Gibbs samples statistics 
(mean, median and credibility intervals). 
'''

import sys
import os
import subprocess
import gzip
import numpy

from utils import *


printScriptHeader()

if len(sys.argv) != 3:

  print("Usage: python calc_rpvg_gibbs_stats.py <rpvg_gibbs_gz_name> <output_file_name>\n")
  sys.exit(1)


read_count_sums = []

gibbs_file = gzip.open(sys.argv[1], "rb")

for line in gibbs_file:

    line_split = line.decode().strip().split("\t")

    if line_split[0] == "Name":

        assert(len(line_split) > 2)
        assert(len(read_count_sums) == 0)      

        for i in range(0, len(line_split) - 2):

            read_count_sums.append(0)

        continue

    assert(len(line_split) == len(read_count_sums) + 2)

    for i in range(0, len(read_count_sums)):

        read_count_sums[i] = read_count_sums[i] + float(line_split[i + 2])

gibbs_file.close()

print(min(read_count_sums))
print(numpy.mean(read_count_sums))
print(numpy.median(read_count_sums))
print(max(read_count_sums))

out_file = open(sys.argv[2], "w")
out_file.write("Name\tClusterID\tMean\tMedian\tMin\tMax\tCILower95\tCIUpper95\tCILower90\tCIUpper90\n")

gibbs_file = gzip.open(sys.argv[1], "rb")

for line in gibbs_file:

    line_split = line.decode().strip().split("\t")

    if line_split[0] == "Name":

        continue

    assert(len(line_split) == len(read_count_sums) + 2)

    samples = []

    for i in range(0, len(read_count_sums)):

        samples.append(float(line_split[i + 2]) / read_count_sums[i] * 1000000)

    samples.sort()

    cil95 = samples[int(round(0.025  * len(samples)))]
    ciu95 = samples[int(round(0.975 * len(samples)))]
    cil90 = samples[int(round(0.05 * len(samples)))]
    ciu90 = samples[int(round(0.95 * len(samples)))]

    out_line = [line_split[0], line_split[1], numpy.mean(samples), numpy.median(samples), min(samples), max(samples), cil95, ciu95, cil90, ciu90]
    out_file.write("\t".join([str(x) for x in out_line]) + "\n")

gibbs_file.close()
out_file.close()

print("\nDone")
