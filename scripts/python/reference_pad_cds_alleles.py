
'''
reference_pad_cds_alleles.py
Add reference seqeunces to the flanks of CDS 
alleles from a single gene.   
'''

import sys
import os
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO

from utils import *


def parse_chromosome(filename, chrom_name):

	sequence = ""

	for record in SeqIO.parse(filename, "fasta"):

		if record.id == chrom_name:

			sequence += str(record.seq)

	return sequence


def parse_gene_coords(filename, gene_name):

	gene_coords = ["", "", -1, -1, -1, -1]

	transcript_file = open(filename, "r")

	for line in transcript_file:

		if line[0] == "#":

			continue

		line_split = line.split("\t")
		attributes_split = line_split[8].split(";")

		cur_gene_name = ""
		tags = []
		transcript_type = ""

		for attribute in attributes_split:

			attribute = attribute.strip()

			if attribute[:9] == "gene_name":

				assert(cur_gene_name == "")
				cur_gene_name = attribute.split('"')[1]
	
			if attribute[:3] == "tag":

				tags.append(attribute.split('"')[1])

			if attribute[:15] == "transcript_type":

				assert(transcript_type == "")
				transcript_type = attribute.split('"')[1]

		assert(cur_gene_name != "")

		if cur_gene_name != gene_name:

			continue

		if gene_coords[0] == "":

			gene_coords[0] = line_split[0]

		assert(gene_coords[0] == line_split[0])

		if gene_coords[1] == "":

			gene_coords[1] = line_split[6]

		assert(gene_coords[1] == line_split[6])

		if line_split[2] == "transcript":

			if gene_coords[2] == -1:

				assert(gene_coords[5] == -1)

				gene_coords[2] = int(line_split[3])
				gene_coords[5] = int(line_split[4])

			else:

				gene_coords[2] = min(gene_coords[2], int(line_split[3]))
				gene_coords[5] = max(gene_coords[5], int(line_split[4]))

		elif line_split[2] == "start_codon" and "basic" in tags and transcript_type == "protein_coding":

			if gene_coords[1] == "-":

				line_split[3] = line_split[4]

			if gene_coords[3] == -1:

				gene_coords[3] = int(line_split[3])

			elif gene_coords[3] != int(line_split[3]):

				print("Warning different start codon:")
				print(gene_coords[3])
				print(int(line_split[3]))

		elif line_split[2] == "stop_codon" and "basic" in tags and transcript_type == "protein_coding":

			if gene_coords[1] == "-":

				line_split[4] = line_split[3]

			if gene_coords[4] == -1:

				gene_coords[4] = int(line_split[4])

			elif gene_coords[4] != int(line_split[4]):

				print("Warning different stop codon:")
				print(gene_coords[4])
				print(int(line_split[4]))

	assert(gene_coords[0] != "")
	assert(gene_coords[1] != "")
	assert(not -1 in gene_coords[2:])

	if gene_coords[1] == "+":

		assert(gene_coords[2] <= gene_coords[3])
		assert(gene_coords[3] < gene_coords[4])
		assert(gene_coords[4] <= gene_coords[5])

	else:

		assert(gene_coords[1] == "-")

		gene_coords[2], gene_coords[5] = gene_coords[5], gene_coords[2]


		assert(gene_coords[2] >= gene_coords[3])
		assert(gene_coords[3] > gene_coords[4])
		assert(gene_coords[4] >= gene_coords[5])

	return gene_coords


printScriptHeader()

if len(sys.argv) != 7:

	print("Usage: python reference_pad_cds_alleles.py <cds_alleles_input_name> <genome_fasta_name> <transcripts_gtf_name> <gene_name> <gene_flank_size> <output_fasta_name>\n")
	sys.exit(1)


gene_coords = parse_gene_coords(sys.argv[3], sys.argv[4])
print(gene_coords)

chrom_seq = parse_chromosome(sys.argv[2], gene_coords[0])
print len(chrom_seq)

cds_file = open(sys.argv[1], "r")
out_file = open(sys.argv[6], "w")

gene_flank_size = int(sys.argv[5])

for line in cds_file:

	line = line.strip()	
	line_split = line.split("\t")

	assert(len(line_split) == 2)

	if line_split[0] == "allele":

		continue

	if gene_coords[1] == "+":

		left_flank = chrom_seq[(gene_coords[2] - gene_flank_size - 1):(gene_coords[3] - 1)]
		right_flank = chrom_seq[gene_coords[4]:(gene_coords[5] + gene_flank_size - 1)]

	else:

		assert(gene_coords[1] == "-")

		left_flank = chrom_seq[gene_coords[3]:(gene_coords[2] + gene_flank_size - 1)]
		right_flank = chrom_seq[(gene_coords[5] - gene_flank_size - 1):(gene_coords[4] - 1)]

		left_flank = Seq(left_flank)
		left_flank = str(left_flank.reverse_complement())

		right_flank = Seq(right_flank)
		right_flank = str(right_flank.reverse_complement())

	out_file.write(">" + line_split[0] + "\n")
	out_file.write(left_flank + line_split[1] + right_flank + "\n")

cds_file.close()
out_file.close()

print("Done")
