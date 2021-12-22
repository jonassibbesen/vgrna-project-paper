set -e

# Set file name prefixes
OUT_PREFIX="${MAPPER}_${REF}_real_r1_${REAL}"

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
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/indexes/${REF}/ . --recursive --no-progress

# Use default HISAT2
if [ "${MAPPER}" = "hisat2" ]; then

	# Map reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t -x ${REF}_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S ${OUT_PREFIX}.sam"

	# Compress alignments
	/usr/bin/time -v bash -c "samtools view -O BAM --threads ${CPU} ${OUT_PREFIX}.sam > ${OUT_PREFIX}.bam"

	# Sort and index alignments
	/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}.bam > ${OUT_PREFIX}_sort.bam; mv ${OUT_PREFIX}_sort.bam ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam"
fi

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/real_r1/${REAL}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.sam" --no-progress
