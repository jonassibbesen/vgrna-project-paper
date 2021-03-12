set -e

# Set file name prefixes
OUT_PREFIX="${MAPPER}_${REF}_real_${REAL}"

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
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/indexes/${REF}/ . --recursive --no-progress

# Use default STAR
if [ "${MAPPER}" = "star" ]; then

	# Map reads
	/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir . --readFilesCommand zcat --readFilesIn reads_1.fq.gz reads_2.fq.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_; mv ${OUT_PREFIX}_Aligned.out.bam ${OUT_PREFIX}.bam"

	# Sort and index alignments
	/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}.bam > ${OUT_PREFIX}_sort.bam; mv ${OUT_PREFIX}_sort.bam ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam"
fi

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/real/${REAL}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.sam" --no-progress
