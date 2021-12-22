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

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.xg" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${REF}/ . --recursive --no-progress

# Get gcsa index name
GCSA=$(ls ${REF}_index*.gcsa) 

# Use faster vg map
if [ "${MAPPER}" = "map_fast" ]; then

	# Map reads
	/usr/bin/time -v bash -c "vg map -t ${CPU} --try-up-to 16 --mate-rescues 32 -x ${REF}.xg -g ${GCSA} -f reads_1.fq.gz -f reads_2.fq.gz > ${OUT_PREFIX}.gam"

# Use default vg mpmap
elif [ "${MAPPER}" = "mpmap" ]; then

	# Map reads
	/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna --report-allelic-mapq -x ${REF}.xg -g ${GCSA} -d ${REF}_index.dist -f reads_1.fq.gz -f reads_2.fq.gz > ${OUT_PREFIX}.gamp"
fi	

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/real_r1/${REAL}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
