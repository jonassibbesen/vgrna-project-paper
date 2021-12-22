
'''
filter_hla_allele_alignments.py
Filters HLA allele alignments that does not align to the
correct gene.  
'''

import sys
import os
import subprocess

import gzip
import json

from Bio import SeqIO

from utils import *

def parse_allele_descriptions(filename):

	allele_descriptions = {}

	for record in SeqIO.parse(filename, "fasta"):

		description_split = record.description.split(" ")

		assert(not description_split[0] in allele_descriptions)
		allele_descriptions[description_split[0]] = description_split[1]

	return allele_descriptions

def parse_transcript_gene_names(filename):

	transcript_gene_names = {}

	transcript_file = open(filename, "r")

	for line in transcript_file:

		if line[0] == "#":

			continue

		line_split = line.split("\t")

		if line_split[2] != "transcript":

			continue	

		attributes_split = line_split[8].split(";")

		gene_name = ""
		transcript_id = ""

		for attribute in attributes_split:

			attribute = attribute.strip()

			if attribute[:9] == "gene_name":

				assert(gene_name == "")
				gene_name = attribute.split('"')[1]

			if attribute[:13] == "transcript_id":

				assert(transcript_id == "")
				transcript_id = attribute.split('"')[1]
	
		assert(gene_name != "")
		assert(transcript_id != "")
	
		assert(not transcript_id in transcript_gene_names)
		transcript_gene_names[transcript_id] = gene_name

	transcript_file.close()
	return transcript_gene_names


printScriptHeader()

if len(sys.argv) != 5:

	print("Usage: python filter_hla_allele_alignments.py <alignment_input_name> <allele_fasta_name> <transcript_gtf_name> <output_name>\n")
	sys.exit(1)


hla_alleles = parse_allele_descriptions(sys.argv[2])
transcript_names = parse_transcript_gene_names(sys.argv[3])

print(len(hla_alleles))
print(len(transcript_names))

align_file = open(sys.argv[1], "r")

filtered_hla_alleles = {}

for line in align_file:

	json_line = json.loads(line)

	hla_gene_name = hla_alleles[json_line["name"]].split("*")[0]

	for refpos in json_line["refpos"]:

		if (not refpos["name"] in transcript_names):

			continue

		gene_name = transcript_names[refpos["name"]]
			
	 	if gene_name[:4] != "HLA-":

	 		continue

	 	if hla_gene_name == gene_name[4:]:

			if json_line["name"] in filtered_hla_alleles:

				if int(json_line["score"]) > filtered_hla_alleles[json_line["name"]][0]:

					filtered_hla_alleles[json_line["name"]] = [int(json_line["score"]), line]

			else:
			
				filtered_hla_alleles[json_line["name"]] = [int(json_line["score"]), line]

			break

align_file.close()

print(len(filtered_hla_alleles))

output_file = open(sys.argv[4], "w")

for name, value in filtered_hla_alleles.items():

	output_file.write(value[1])

output_file.close()

print("Done")
