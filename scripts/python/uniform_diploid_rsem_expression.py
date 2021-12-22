
'''
uniform_diploid_rsem_expression.py
Creates new expression file where the TPM values are uniform (1) across 
haplotypes. The expression value (TPM) of transcripts which does not 
have a version on both haplotypes are also set to 0.
'''

import sys
import os
import subprocess

from utils import *

def write_uniform_exp(hap1_line_split, hap2_line_split, exp_out_file):

	assert(len(hap1_line_split) != 0 or len(hap2_line_split) != 0)

	uniform_tpm = 0

	if len(hap1_line_split) != 0 and len(hap2_line_split) != 0:

		uniform_tpm = 1

	if len(hap1_line_split) != 0:

		hap1_line_split[4:] = ["0"] * 4
		hap1_line_split[5] = str(uniform_tpm)

		exp_out_file.write("\t".join(hap1_line_split) + "\n")

	if len(hap2_line_split) != 0:
	
		hap2_line_split[4:] = ["0"] * 4
		hap2_line_split[5] = str(uniform_tpm)

		exp_out_file.write("\t".join(hap2_line_split) + "\n")

printScriptHeader()

if len(sys.argv) != 3:

	print("Usage: python uniform_diploid_rsem_expression.py <input_name> <output_name>\n")
	sys.exit(1)

exp_in_file = open(sys.argv[1], "r")
exp_out_file = open(sys.argv[2], "w")

hap1_line_split = []
hap2_line_split = []

prev_transcript_id = ""

for line in exp_in_file:

	line_split = line.strip().split("\t")
	
	if line_split[0] == "transcript_id":

		exp_out_file.write(line)
		continue

	# Assumes haplotype-specific transcripts from the same transcript are 
	# after each other.
	assert(line_split[0][-2] == "_" or hap_id == "_")

	hap_id = line_split[0][-1]
	cur_transcript_id = line_split[0][:-2]

	if cur_transcript_id != prev_transcript_id and prev_transcript_id != "":

		write_uniform_exp(hap1_line_split, hap2_line_split, exp_out_file)

		hap1_line_split = []
		hap2_line_split = []

	if hap_id == "1":

		assert(len(hap1_line_split) == 0)
		hap1_line_split = line_split

	else:

		assert(hap_id == "2")
		assert(len(hap2_line_split) == 0)
		hap2_line_split = line_split

	prev_transcript_id = cur_transcript_id

write_uniform_exp(hap1_line_split, hap2_line_split, exp_out_file)

exp_in_file.close()
exp_out_file.close()

print("Done")
