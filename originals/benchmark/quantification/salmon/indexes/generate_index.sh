set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"
OUT_PREFIX="salmon_index_${INDEX}"

# Download genome
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa.fai . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.fa.gz" --exclude "*MT*" --exclude "*SCA*" --no-progress

if echo "${INDEX}" | grep -q "decoy"; then 

	# Combine transcripts
	/usr/bin/time -v bash -c "cat */*fa.gz <(gzip -c ${GENOME_PREFIX}.fa) > ${REF}.fa.gz; zcat ${REF}.fa.gz | grep '"'>'"' | wc -l"

	# Generate decoy list
	/usr/bin/time -v bash -c "grep '"'>'"' ${GENOME_PREFIX}.fa | cut -d ' ' -f 1 | sed -e 's/>//g' > decoys.txt; wc -l decoys.txt"

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -d decoys.txt -i ${OUT_PREFIX} -t ${REF}.fa.gz"

else 

	# Combine transcripts
	/usr/bin/time -v bash -c "cat */*fa.gz > ${REF}.fa.gz; zcat ${REF}.fa.gz | grep '"'>'"' | wc -l"

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -i ${OUT_PREFIX} -t ${REF}.fa.gz"
fi

# Upload index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/salmon/indexes/${INDEX}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
