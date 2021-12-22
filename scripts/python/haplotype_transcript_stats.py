
'''
haplotype_transcript_stats.py
Calculate stats on haplotype-specific transcripts shared
between samples, populations and super populations.
'''

import sys
import gzip

from utils import *

printScriptHeader()

if len(sys.argv) != 3:

	print("Usage: python haplotype_transcript_stats.py <1kg_pop_info> <hap_txp_info> > stats.txt\n")
	sys.exit(1)


pop_file = open(sys.argv[1], "r") 
pop_data = {}

for line in pop_file:

	line_split = line.split("\t")
	assert(len(line_split) >= 4)

	if line_split[0] != "sample":

		pop_data[line_split[0]] = (line_split[1], line_split[2])

pop_file.close()

sys.stderr.write(str(len(pop_data)) + "\n")


hap_file = open(sys.argv[2], "r") 
hap_data = []
hap_index = {}

for line in hap_file:

	line_split = line.split("\t")
	assert(len(line_split) == 5)

	if line_split[0] != "Name":

		hap_data.append([x.split("_")[2] for x in line_split[4].split(",")])

		for hap in hap_data[-1]:

			if hap in hap_index:

				hap_index[hap].append(len(hap_data) - 1)

			else:

				hap_index[hap] = [len(hap_data) - 1]


hap_file.close()

sys.stderr.write(str(len(hap_data)) + "\n")
sys.stderr.write(str(len(hap_index)) + "\n")


print("name\tpop\tspop\ttotal\tnum_sam\tnum_pop\tnum_spop")

for pop_key, pop_value in pop_data.items():

	num_haps = 0
	num_sam = 0
	num_pop = 0
	num_spop = 0

	for hap_idx in hap_index[pop_key]:

		num_haps += 1

		has_sam = False
		has_pop = False
		has_spop = False

		for hap in hap_data[hap_idx]:

			pop = pop_data[hap]

			if not has_sam and hap != pop_key:

				num_sam += 1
				has_sam = True

			if not has_pop and pop[0] != pop_value[0]:

				num_pop += 1
				has_pop = True

			if not has_spop and pop[1] != pop_value[1]:

				num_spop += 1
				has_spop = True

			if has_sam and has_pop and has_spop:

				break


	print(pop_key + "\t" + pop_value[0] + "\t" + pop_value[1] + "\t" + str(num_haps) + "\t" + str(num_sam) + "\t" + str(num_pop) + "\t" + str(num_spop))
	sys.stdout.flush()


