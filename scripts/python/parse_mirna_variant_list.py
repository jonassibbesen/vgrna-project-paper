
'''
parse_mirna_variant_list.py
Parses and filters miRNA variant list from Fehlmann et at.
paper. Output type is vcf.  
'''

import sys
import os
import subprocess
import random

from Bio import SeqIO

from utils import *

def parse_annotation(anno_filename):

	anno_file = open(anno_filename, "r")

	annotation = {}

	for line in anno_file:

		if line[0] == "#":

			continue

		line_split = line.strip().split("\t")

		if line_split[2] != "miRNA":

			continue

		att_split = line_split[8].split(";")
		assert(att_split[2][:4] == "Name")

		name = att_split[2][5:]

		if name in annotation:

			annotation[name].append((line_split[0], int(line_split[3]), int(line_split[4]), line_split[6]))

		else:

			annotation[name] = [(line_split[0], int(line_split[3]), int(line_split[4]), line_split[6])]

	anno_file.close()

	return annotation

def parse_genome(genome_filename):

	genome = {}

	for record in SeqIO.parse(genome_filename, "fasta"):

		if record.id in genome:

			genome[record.id] += str(record.seq)

		else:

			genome[record.id] = str(record.seq)

	return genome

printScriptHeader()

if len(sys.argv) != 6:

	print("Usage: python parse_mirna_variant_list.py <mirna_variants_csv> <mirna_annotation_gff> <genome_fasta> <min_relative_count> <output_name>\n")
	sys.exit(1)

variant_in_file = open(sys.argv[1], "r")

mirnas = parse_annotation(sys.argv[2])
print(len(mirnas))

genome = parse_genome(sys.argv[3])
print(len(genome))

min_relative_count = float(sys.argv[4])
assert(min_relative_count >= 0)

mirna_variants = {}

for line in variant_in_file:

	line_split = line.strip().split(";")
	assert(len(line_split) == 6)
	
	if line_split[0] == "ptecutsot":

		continue

	name = line_split[1]

	if name in mirna_variants:

		if line_split[2] in mirna_variants[name]:

			mirna_variants[name][line_split[2]] == max(mirna_variants[name][line_split[2]], int(line_split[5]))

		else:

			mirna_variants[name][line_split[2]] = int(line_split[5])

	else:

		mirna_variants[name] = {}
		mirna_variants[name][line_split[2]] = int(line_split[5])

variant_in_file.close()

print(len(mirna_variants))

bla = 0

for key, value in mirna_variants.items():

	if key in mirnas:

		if len(value) < 20:

			mir = mirnas[key][0]

			if mir[3] == "+":

				print(key)
				print(value)
				print(mirnas[key])


				print(genome[mir[0]][(mir[1] - 1):mir[2]])

				for key2, value2 in value.items():

					if len(key2) == 0:

						continue

					print(key2 + "\t" + genome[mir[0]][(mir[1] - 1) + int(key2.split(":")[0])])
		

	else:

		print("Warning: " + key + " not in miRNA annotation")



vcf_out_file = open(sys.argv[5], "w")



vcf_out_file.close()

print("Done")

