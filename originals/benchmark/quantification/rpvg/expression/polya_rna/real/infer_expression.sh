set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${GRAPH}_real_${REAL}"
OUT_PREFIX="${QUANTER}_${MAPPER}_${TRANSCRIPTS}_real_${REAL}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/real/${REAL}/${MAPPER}/${GRAPH}/ . --recursive --exclude "*" --include "*.gam" --include "*.gamp" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt.ri . --no-progress

# Download transcript info
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Concatenate transcript info
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep -v ^Name > transcript_info.txt; wc -l transcript_info.txt"

if [ "${MAPPER}" = "map" ]; then 

	if [ "${QUANTER}" = "rpvg_strand" ]; then 

		# Infer expression
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"
	fi

elif [ "${MAPPER}" = "mpmap" ] || [ "${MAPPER}" = "mpmap_nosplice" ]; then

	if [ "${QUANTER}" = "rpvg" ]; then 

		# Infer expression
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_strand" ]; then 

		# Infer expression
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_strand_gam" ]; then

		# Convert to gamp to gam
		/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"

		# Infer expression        
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"
	fi
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rpvg/expression/polya_rna/real/${REAL}/${QUANTER}/${TRANSCRIPTS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
