set -e

# Set file name prefixes
OUT_PREFIX="${QUANTER}_${REF}_real_r2_${REAL}"

# Reads are in local aws bucket
if (echo "${READS_1}" | grep -q "jsibbesen"); then 

	# Download reads
	aws s3 cp ${READS_1} reads_1.fq.gz --no-progress
	aws s3 cp ${READS_2} reads_2.fq.gz --no-progress

else

	# Download reads
	aws s3 cp ${READS_1} reads_1.fq.gz --no-sign-request --no-progress
	aws s3 cp ${READS_2} reads_2.fq.gz --no-sign-request --no-progress
fi

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/indexes/${REF}/ . --recursive --no-progress

# Use stranded RSEM
if [ "${QUANTER}" = "rsem_strand" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --no-bam-output --strandedness reverse --paired-end reads_1.fq.gz reads_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"

# Use stranded RSEM with more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_strand_k1k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 1000 --seed ${SEED} --no-bam-output --strandedness reverse --paired-end reads_1.fq.gz reads_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"

# Use stranded RSEM with many more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_strand_k2k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 2000 --seed ${SEED} --no-bam-output --strandedness reverse --paired-end reads_1.fq.gz reads_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.genes.results; gzip ${OUT_PREFIX}.isoforms.results"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/expression/polya_rna/real_r2/${REAL}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
