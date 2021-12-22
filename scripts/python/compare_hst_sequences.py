
'''
compare_hst_sequences.py
Compares haplotype-specific transcript (HST) sequences. Only the same
transcripts are compared. Writes a tsv file with the names of the 
seqeunces that are identical.  
'''

import sys
import os
import subprocess
import random

from Bio import SeqIO

from utils import *

def parse_transcripts(filename):

	transcripts = {}

	for record in SeqIO.parse(filename, "fasta"):

		transcript_id = record.id.split("_")[0]

		if transcript_id in transcripts:

			transcripts[transcript_id].append((record.id, str(record.seq)))

		else:

			transcripts[transcript_id] = [(record.id, str(record.seq))]

	return transcripts


printScriptHeader()

if len(sys.argv) != 4:

	print("Usage: python compare_hst_sequences.py <input_name_1> <input_name_2> <output_name>\n")
	sys.exit(1)


hts_seqs_1 = parse_transcripts(sys.argv[1])
print("Parsed " + str(len(hts_seqs_1)) + " transcripts")

hts_seqs_2 = parse_transcripts(sys.argv[2])
print("Parsed " + str(len(hts_seqs_2)) + " transcripts")

tsv_out_file = open(sys.argv[3], "w")
tsv_out_file.write("Name1\tName2\n")

for name_1, seqs_1 in hts_seqs_1.items():

	if name_1 in hts_seqs_2:

		for hts_2 in hts_seqs_2[name_1]:

			for hts_1 in seqs_1:

				if hts_1[1] == hts_2[1]:

					tsv_out_file.write(hts_1[0] + "\t" + hts_2[0] + "\n")

tsv_out_file.close()

print("Done")
