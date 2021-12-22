
'''
calc_vg_sim_fragment_length_stats.py
Calculates fragment length mean and standard deviation 
from vg sim path position output.
'''

import sys
import gzip
import numpy

from utils import *

printScriptHeader()

if len(sys.argv) != 3:

	print("Usage: python calc_vg_sim_fragment_length_stats.py <input_name> <read_length>\n")
	sys.exit(1)


pos_file = gzip.open(sys.argv[1], "rb")
read_length = int(sys.argv[2])

frag_lengths = {}
unpaired_mates = {}

num_reads = 0

for line in pos_file:

	line_split = line.split("\t")
	assert(len(line_split) == 4)

	if line_split[0] != "read":

		num_reads += 1

		read_split = line_split[0].split("_")
		assert(len(read_split) == 5)

		if read_split[3] in unpaired_mates:

			cur_pos1 = unpaired_mates.pop(read_split[3])
			cur_length = abs(cur_pos1 - int(line_split[2])) + read_length

			if cur_length in frag_lengths:

				frag_lengths[cur_length] += 1

			else:

				frag_lengths[cur_length] = 1

		else:

			unpaired_mates[read_split[3]] = int(line_split[2])

	if num_reads % 10000000 == 0:

		print num_reads

pos_file.close()

print num_reads
print(frag_lengths)

total_count = 0.0;
sum_count = 0.0; 

for key, value in frag_lengths.items():

    total_count += value
    sum_count += (key * value)

mean = (sum_count / total_count);

sum_var = 0.0; 

for key, value in frag_lengths.items():

	sum_var += (((float(key) - mean) ** 2) * value);
  

sd = ((sum_var / (total_count - 1)) ** 0.5);

print("mean: ", mean)
print("sd: ", sd)



