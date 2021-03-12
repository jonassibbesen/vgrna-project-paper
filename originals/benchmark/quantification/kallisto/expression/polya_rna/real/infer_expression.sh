set -e

# Set file name prefixes
OUT_PREFIX="${QUANTER}_${REF}_real_${REAL}"

# Reads are in aws s3 bucket
if (echo "${READS_1}" | grep -q "s3:"); then 

	# Reads are in local bucket
	if (echo "${READS_1}" | grep -q "jsibbesen"); then 

		# Download reads
		aws s3 cp ${READS_1} reads_1.fq.gz --no-progress
		aws s3 cp ${READS_2} reads_2.fq.gz --no-progress

	else

		# Download reads
		aws s3 cp ${READS_1} reads_1.fq.gz --no-sign-request --no-progress
		aws s3 cp ${READS_2} reads_2.fq.gz --no-sign-request --no-progress
	fi

else

	# Download reads
	wget --no-check-certificate --no-verbose -O reads_1.fq.gz ${READS_1}
	wget --no-check-certificate --no-verbose -O reads_2.fq.gz ${READS_2}
fi

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/kallisto/indexes/${REF}/ . --recursive --no-progress

# Use Kallisto
if [ "${QUANTER}" = "kallisto" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} -i kallisto_index_${REF}.idx -o ${OUT_PREFIX} reads_1.fq.gz reads_2.fq.gz"

# Use stranded Kallisto
elif [ "${QUANTER}" = "kallisto_strand" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} --rf-stranded -i kallisto_index_${REF}.idx -o ${OUT_PREFIX} reads_1.fq.gz reads_2.fq.gz"

# Use Salmon with bias correction
elif [ "${QUANTER}" = "kallisto_strand_bias" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} --rf-stranded --bias -i kallisto_index_${REF}.idx -o ${OUT_PREFIX} reads_1.fq.gz reads_2.fq.gz"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/abundance.tsv"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/kallisto/expression/polya_rna/real/${REAL}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
