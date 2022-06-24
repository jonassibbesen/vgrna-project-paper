
'''
convert_salmon_bootstraps.py 
Convert Salmon bootstrap values to mean 
expression values. Inspired by 
https://github.com/COMBINE-lab/salmon/blob/master/scripts/ConvertBootstrapsToTSV.py
'''

import sys
import os
import subprocess
import gzip
import struct

from utils import *


def parse_hst_names(filename):

    hst_names = []
    
    names_file = gzip.open(filename, "rb")

    for line in names_file:

        assert(len(hst_names) == 0)
        hst_names = line.decode().strip().split("\t")

    names_file.close()

    return hst_names

def parse_hst_lengths(filename):

    hst_lengths = {}
    
    lengths_file = gzip.open(filename, "rb")

    for line in lengths_file:

        line_split = line.decode().split("\t")
        assert(len(line_split) == 5)

        if line_split[0] == "Name":

            continue

        assert(not line_split[0] in hst_lengths)
        hst_lengths[line_split[0]] = (int(line_split[1]), float(line_split[2]))

    lengths_file.close()

    return hst_lengths


printScriptHeader()

if len(sys.argv) != 5:

    print("Usage: python convert_salmon_bootstraps.py <bootstraps_gz_name> <names_tsv_gz_name> <em_quant_gz_name> <output_file_name>\n")
    sys.exit(1)


hst_names = parse_hst_names(sys.argv[2])
print(len(hst_names))

hst_lengths = parse_hst_lengths(sys.argv[3])
print(len(hst_lengths))

boot_struct = struct.Struct('@' + 'd' * len(hst_names))

sum_boot_values = [0] * len(hst_names)
num_zero_boot_values = [0] * len(hst_names)

num_boot_samples = 0

with gzip.open(sys.argv[1], "rb") as boot_file:

    while True:

        try:

            boot_values = boot_struct.unpack_from(boot_file.read(boot_struct.size))    
            assert(len(sum_boot_values) == len(boot_values))

            for i in range(len(sum_boot_values)):

                sum_boot_values[i] += float(boot_values[i])

                if float(boot_values[i]) < 10**-4:

                    num_zero_boot_values[i] += 1
        
            num_boot_samples += 1

        except:

            break

boot_file.close()

print(num_boot_samples)

mean_boot_count = [0] * len(hst_names)
mean_boot_tpm = [0] * len(hst_names)

for i in range(len(sum_boot_values)):

    mean_boot_count[i] = sum_boot_values[i] / num_boot_samples
    mean_boot_tpm[i] = mean_boot_count[i] / hst_lengths[hst_names[i]][1]

assert(len(hst_names) == len(mean_boot_count))
assert(len(hst_names) == len(mean_boot_tpm))
assert(len(hst_names) == len(num_zero_boot_values))

sum_mean_boot_tpm = sum(mean_boot_tpm)

out_file = open(sys.argv[4], "w")
out_file.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\tFracNonZero\n")

for i in range(len(hst_names)):

    out_str = hst_names[i] + "\t" + str(hst_lengths[hst_names[i]][0]) + "\t" + str(hst_lengths[hst_names[i]][1]) + "\t" + str(10**6 * mean_boot_tpm[i] / sum_mean_boot_tpm) + "\t" + str(mean_boot_count[i]) + "\t" + str(1 - num_zero_boot_values[i] / num_boot_samples) + "\n"

    out_file.write(out_str)

out_file.close()

print("Done")
