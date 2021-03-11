set -e

# Set genome prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set transcripts prefix
TRANSCRIPTS_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="salmon_index_1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

if echo "${INDEX}" | grep -q "decoy"; then 

	# Generate decoy list
	/usr/bin/time -v bash -c "grep '"'>'"' ${GENOME_PREFIX}.fa | cut -d ' ' -f 1 | sed -e 's/>//g' > decoys.txt; wc -l decoys.txt"

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -d decoys.txt -i ${OUT_PREFIX} -t ${TRANSCRIPTS_PREFIX}.fa.gz; rm "

else 

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -i ${OUT_PREFIX} -t ${TRANSCRIPTS_PREFIX}.fa.gz"
fi
