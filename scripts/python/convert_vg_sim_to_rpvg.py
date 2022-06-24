
'''
convert_vg_sim_to_rpvg.py 
Convert vg sim path count output to rpvg 
expression output.  
'''

import sys
import os
import subprocess

import pickle
import gzip

from Bio.Seq import Seq
from Bio import SeqIO

from utils import *


def parse_path_counts(filename, is_paired):

	path_counts = {}
	
	vg_sim_file = gzip.open(filename, "rb")

	for line in vg_sim_file:

		line_split = line.strip().split("\t")
		assert(len(line_split) == 4)

		if line_split[0] == "read":

			continue

		if line_split[1] in path_counts: 

			if is_paired:

				path_counts[line_split[1]] += 0.5

			else:

				path_counts[line_split[1]] += 1				

		else:

			if is_paired:

				path_counts[line_split[1]] = 0.5

			else:

				path_counts[line_split[1]] = 1


	vg_sim_file.close()

	return path_counts

def parse_isoforms_lengths(filename):

	isoform_lengths = {}
	
	isoforms_file = open(filename, "rb")

	for line in isoforms_file:

		line_split = line.strip().split("\t")
		assert(len(line_split) == 8)

		if line_split[0] == "transcript_id":

			continue

		assert(not line_split[0] in isoform_lengths)
		isoform_lengths[line_split[0]] = [int(line_split[2]), float(line_split[3])]

	isoforms_file.close()

	return isoform_lengths


printScriptHeader()

if len(sys.argv) != 5:

	print("Usage: python convert_vg_sim_to_rpvg.py <vg_sim_gz_name> <is_paired (Y|N)> <isoform_length_name> <output_file_name>\n")
	sys.exit(1)


assert(sys.argv[2] == "Y" or sys.argv[2] == "N")

path_counts = parse_path_counts(sys.argv[1], sys.argv[2] == "Y")
print(len(path_counts))

isoform_lengths = parse_isoforms_lengths(sys.argv[3])

total_transcript_count = 0

for path, count in path_counts.items():

	total_transcript_count += (count / isoform_lengths[path][1])

print(total_transcript_count)

out_file = open(sys.argv[4], "w")
out_file.write("Name\tClusterID\tLength\tEffectiveLength\tHaplotypeProbability\tReadCount\tTPM\n")	

for path, count in path_counts.items():

	length = isoform_lengths[path]
	tpm = (count / length[1]) * 10**6 / total_transcript_count

	out_file.write(path + "\t0\t" + str(length[0]) + "\t" + str(length[1]) + "\t1\t" + str(count) + "\t" + str(tpm) + "\n")

out_file.close()

print("Done")
