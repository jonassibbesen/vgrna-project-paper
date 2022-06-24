set -e

# Set file name prefixes
OUT_PREFIX="${QUANTER}_${REF}_real_r2_${REAL}"

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
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/salmon/indexes/${REF}/ . --recursive --no-progress

# Use default Salmon
if [ "${QUANTER}" = "salmon" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A -i salmon_index_${REF} -o ${OUT_PREFIX} -1 reads_1.fq.gz -2 reads_2.fq.gz"

# Use Salmon with bias correction
elif [ "${QUANTER}" = "salmon_bias" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --seqBias --gcBias -i salmon_index_${REF} -o ${OUT_PREFIX} -1 reads_1.fq.gz -2 reads_2.fq.gz"

# Use Salmon with traditional EM algorithm
elif [ "${QUANTER}" = "salmon_em" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --useEM -i salmon_index_${REF} -o ${OUT_PREFIX} -1 reads_1.fq.gz -2 reads_2.fq.gz"

# Use Salmon with traditional EM algorithm and bias correction
elif [ "${QUANTER}" = "salmon_em_bias" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --useEM --seqBias --gcBias -i salmon_index_${REF} -o ${OUT_PREFIX} -1 reads_1.fq.gz -2 reads_2.fq.gz"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/quant.sf"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/salmon/expression/polya_rna/real_r2/${REAL}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
