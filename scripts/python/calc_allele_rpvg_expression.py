
'''
calc_allele_rpvg_expression.py 
Calculates allele expression from rpvg expression
estimates. Note that the input variant (vcf) file 
need to annotated with transcript names 
(INFO:TRANSCIPTS tag).  
'''

import sys
import os
import subprocess

import pickle
import gzip

from Bio.Seq import Seq
from Bio import SeqIO

from utils import *


def parse_hst_info(filename):

	hst_info = {}
	
	hst_file = gzip.open(filename, "rb")

	for line in hst_file:

		line_split = line.split("\t")
		assert(len(line_split) == 5)

		if line_split[0] == "Name":

			continue

		if not line_split[2] in hst_info:

			hst_info[line_split[2]] = []

		hst_info[line_split[2]].append((line_split[0], [line_split[4].split(",")[0].split("_")[i] for i in [2,4]]))

	hst_file.close()

	return hst_info


def parse_genome(filename):

	genome = {}

	for record in SeqIO.parse(filename, "fasta"):

		if not record.id in genome: 

			genome[record.id] = ""

		genome[record.id] += str(record.seq)

	return genome


def parse_rpvg_haplotypes(filename):

	rpvg_haps = {}
	
	rpvg_file = gzip.open(filename, "rb")

	for line in rpvg_file:

		line_split = line.split("\t")
		assert(len(line_split) >= 4)

		if line_split[0] == "Name1" or line_split[0] == "Name_1":

			continue

		if not line_split[0] in rpvg_haps:

			rpvg_haps[line_split[0]] = {}

		assert(not line_split[1] in rpvg_haps[line_split[0]])
		rpvg_haps[line_split[0]][line_split[1]] = float(line_split[3])

		if line_split[0] != line_split[1]:
	
			if not line_split[1] in rpvg_haps:

				rpvg_haps[line_split[1]] = {}

			assert(not line_split[0] in rpvg_haps[line_split[1]])
			rpvg_haps[line_split[1]][line_split[0]] = float(line_split[3])

	rpvg_file.close()

	return rpvg_haps


def parse_rpvg_expression(filename, frag_seq_length):

	rpvg_exp = {}
	
	rpvg_file = gzip.open(filename, "rb")

	for line in rpvg_file:

		line_split = line.strip().split("\t")
		assert(len(line_split) == 7)

		if line_split[0] == "Name":

			continue

		assert(not line_split[0] in rpvg_exp)

		if line_split[0] != "Unknown":

			assert((float(line_split[5]) > 0) == (float(line_split[6]) > 0));

			if float(line_split[5]) > 0:

				rpvg_exp[line_split[0]] = [float(line_split[5]) * frag_seq_length / float(line_split[2]), float(line_split[5]), float(line_split[6])]

	rpvg_file.close()

	return rpvg_exp


printScriptHeader()

if len(sys.argv) != 8:

	print("Usage: python calc_allele_rpvg_expression.py <variant_vcf_gz_name> <hst_input_gz_name> <genome_fasta_file> <rpvg_haplotypes_gz_name> <rpvg_expression_gz_name> <frag_seq_length> <output_file_name>\n")
	sys.exit(1)


hst_info = parse_hst_info(sys.argv[2])
print(len(hst_info))

genome = parse_genome(sys.argv[3])
print(len(genome))

if sys.argv[4] == "null":

	rpvg_haps = {}

else:

	rpvg_haps = parse_rpvg_haplotypes(sys.argv[4])

print(len(rpvg_haps))

rpvg_exp = parse_rpvg_expression(sys.argv[5], int(sys.argv[6]))
print(len(rpvg_exp))

out_file = open(sys.argv[7], "w")
out_file.write("Chrom\tPosition\tAlleleSeq\tAlleleNum\tAlleleType\tAlleleLength\tHomopolymerLength\tNumTandemRepeats\tProbability\tBaseReadCount\tTranscriptReadCount\tTPM\n")	

variant_file = gzip.open(sys.argv[1], "rb")

sample_names = {}

for line in variant_file:

	line_split = line.split("\t")
	line_split[-1] = line_split[-1].strip()

	if line_split[0] == "#CHROM":

		assert(len(line_split) >= 10)

		for i in range(9, len(line_split)):

			sample_names[line_split[i]] = i

		continue

	if line_split[0][0] == "#":

		continue

	assert(len(line_split) >= 10)

	alt_alleles = line_split[4].split(",")

	allele_prob = [0.0 for x in range(1 + len(alt_alleles))]	
	allele_exp = [[0.0, 0.0, 0.0] for x in range(1 + len(alt_alleles))]

	transcripts = [x.split("=")[1].split(",") for x in line_split[7].split(";") if x.split("=")[0] == "TRANSCIPTS"]
	assert(len(transcripts) == 1)

	for transcript in transcripts[0]:

		added_diplotypes = [{} for x in range(1 + len(alt_alleles))]

		if transcript in hst_info:

			for hst in hst_info[transcript]:

				gt = line_split[sample_names[hst[1][0]]]
				assert(not "/" in gt)

				allele = gt.split("|")[int(hst[1][1])]

				if allele != ".":

					if hst[0] in rpvg_haps:

						for key, value in rpvg_haps[hst[0]].items():

							if not (hst[0], key) in added_diplotypes[int(allele)] and not (key, hst[0]) in added_diplotypes[int(allele)]:

								added_diplotypes[int(allele)][(hst[0], key)] = ""
								allele_prob[int(allele)] += value

					if hst[0] in rpvg_exp:

						allele_exp[int(allele)][0] += rpvg_exp[hst[0]][0]
						allele_exp[int(allele)][1] += rpvg_exp[hst[0]][1]
						allele_exp[int(allele)][2] += rpvg_exp[hst[0]][2]

	if len(transcripts[0]) > 0:

		for i in range(len(allele_prob)):

			allele_prob[i] = allele_prob[i] / len(transcripts[0])

	max_hp_length = calcMaxHomopolymerLength(genome[line_split[0]], int(line_split[1]) - 1)
	out_file.write(line_split[0] + "\t" + line_split[1] + "\t" + line_split[3] + "\t0\tRef\t0\t" + str(max_hp_length) + "\t0\t" + str(allele_prob[0]) + "\t" + str(allele_exp[0][0]) + "\t" + str(allele_exp[0][1]) + "\t" + str(allele_exp[0][2]) + "\n")

	for i in range(len(alt_alleles)):

		allele_type_length = getAlleleTypeLength(line_split[3], alt_alleles[i])
		assert(allele_type_length[0] != "Ref")

		max_num_tr = calcMaxNumTandemRepeats(genome[line_split[0]], int(line_split[1]) - 1, line_split[3], alt_alleles[i])

		out_file.write(line_split[0] + "\t" + line_split[1] + "\t" + alt_alleles[i] + "\t" + str(i + 1) + "\t" + allele_type_length[0] + "\t" + str(allele_type_length[1]) + "\t" + str(max_hp_length) + "\t" + str(max_num_tr) + "\t" + str(allele_prob[i + 1]) + "\t" + str(allele_exp[i + 1][0]) + "\t" + str(allele_exp[i + 1][1]) + "\t" + str(allele_exp[i + 1][2]) + "\n")

variant_file.close()
out_file.close()

print("Done")
