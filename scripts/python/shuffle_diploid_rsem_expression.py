
'''
shuffle_diploid_rsem_expression.py
Creates new expression file where the TPM values are randomly shuffled
between transcripts. Haplotype-specific transcripts are shuffled seperately  
from transcripts with the same expression value (TPM) between haplotypes.
'''

import sys
import os
import subprocess
import random

from utils import *

class ExpressionValues:

	equal = []
	diff = [] 
	hap1 = []
	hap2 = []	

def add_exp_values(hap1_line_split, hap2_line_split, exp_values):

	assert(len(hap1_line_split) != 0 or len(hap2_line_split) != 0)

	if len(hap1_line_split) != 0 and len(hap2_line_split) != 0:

		hap1_tpm = hap1_line_split[5]
		hap2_tpm = hap2_line_split[5]

		if (round(float(hap1_tpm), 2) == round(float(hap2_tpm), 2)):

			exp_values.equal.append((hap1_tpm, hap1_tpm))

		else:

			exp_values.diff.append((hap1_tpm, hap2_tpm))

	elif len(hap1_line_split) != 0:

		exp_values.hap1.append(hap1_line_split[5])

	else: 

		assert(len(hap2_line_split) != 0)
		exp_values.hap2.append(hap2_line_split[5])

def write_shuffled_exp(hap1_line_split, hap2_line_split, exp_values, exp_out_file):

	assert(len(hap1_line_split) != 0 or len(hap2_line_split) != 0)

	if len(hap1_line_split) != 0 and len(hap2_line_split) != 0:

		if (round(float(hap1_line_split[5]), 2) == round(float(hap2_line_split[5]), 2)):

			new_tpms = exp_values.equal.pop()

		else:

			new_tpms = exp_values.diff.pop()

		hap1_line_split[5] = new_tpms[0]
		hap2_line_split[5] = new_tpms[1]

	elif len(hap1_line_split) != 0:

		hap1_line_split[5] = exp_values.hap1.pop()

	else: 

		assert(len(hap2_line_split) != 0)
		hap2_line_split[5] = exp_values.hap2.pop()

	if len(hap1_line_split) != 0:

		hap1_line_split[4] = "0"
		hap1_line_split[6:] = ["0"] * 2

		exp_out_file.write("\t".join(hap1_line_split) + "\n")

	if len(hap2_line_split) != 0:
		
		hap2_line_split[4] = "0"
		hap2_line_split[6:] = ["0"] * 2

		exp_out_file.write("\t".join(hap2_line_split) + "\n")

printScriptHeader()

if len(sys.argv) != 4:

	print("Usage: python shuffle_diploid_rsem_expression.py <input_name> <output_name> <seed>\n")
	sys.exit(1)

exp_in_file = open(sys.argv[1], "r")
exp_out_file = open(sys.argv[2], "w")

hap1_line_split = []
hap2_line_split = []

prev_transcript_id = ""

exp_values = ExpressionValues() 

for line in exp_in_file:

	line_split = line.strip().split("\t")
	
	if line_split[0] == "transcript_id":

		continue

	# Assumes haplotype-specific transcripts from the same transcript are 
	# after each other.
	assert(line_split[0][-2] == "_" or hap_id == "_")

	hap_id = line_split[0][-1]
	cur_transcript_id = line_split[0][:-2]

	if cur_transcript_id != prev_transcript_id and prev_transcript_id != "":

		add_exp_values(hap1_line_split, hap2_line_split, exp_values)

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

add_exp_values(hap1_line_split, hap2_line_split, exp_values)

random.seed(int(sys.argv[3]))

random.shuffle(exp_values.equal)
random.shuffle(exp_values.diff)
random.shuffle(exp_values.hap1)
random.shuffle(exp_values.hap2)

exp_in_file.seek(0)
prev_transcript_id = ""

hap1_line_split = []
hap2_line_split = []

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

		write_shuffled_exp(hap1_line_split, hap2_line_split, exp_values, exp_out_file)

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

write_shuffled_exp(hap1_line_split, hap2_line_split, exp_values, exp_out_file)

exp_in_file.close()
exp_out_file.close()

print("Done")
