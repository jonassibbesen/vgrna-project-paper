
'''
create_encode_read_data_lists.py
Extract ENCODE meta data and create reads data lists. 
'''

import requests
import json

from utils import *

if len(sys.argv) != 2:

	print("Usage: python create_encode_read_data_lists.py <experiment>\n")
	sys.exit(1)


headers = {"accept": "application/json"}
response = requests.get('https://www.encodeproject.org/experiments/' + str(sys.argv[1]), headers=headers)

experiment = response.json()

assert(experiment["assay_term_name"] == "polyA plus RNA-seq")

print(sys.argv[1])
print(experiment["biosample_ontology"]["term_name"])

reads = {}

for exp_files in experiment["files"]:

	if exp_files["file_format"] != "fastq":

		continue

	assert(exp_files["run_type"] == "paired-ended")

	cur_accession = exp_files["accession"]
	mate_accession = exp_files["paired_with"].split("/")[2]

	if exp_files["paired_end"] == "1":

		if cur_accession in reads:

			assert(reads[cur_accession][1][1] == cur_accession)
			reads[cur_accession][0] = [exp_files["s3_uri"], mate_accession]

		else:

			reads[cur_accession] = [[exp_files["s3_uri"], mate_accession], ["",""]]

	else:

		assert(exp_files["paired_end"] == "2")

		if mate_accession in reads:

			assert(reads[mate_accession][0][1] == cur_accession)
			reads[mate_accession][1] = [exp_files["s3_uri"], mate_accession]

		else:

			reads[mate_accession] = [[exp_files["s3_uri"], mate_accession], ["",""]]

reads_1 = []
reads_2 = []

for key, value in reads.items():

	assert(len(value) == 2)

	assert(len(value[0][0]) != "")
	assert(len(value[1][0]) != "")

	reads_1.append(value[0][0])
	reads_2.append(value[1][0])

print(",".join(reads_1))
print(",".join(reads_2))

