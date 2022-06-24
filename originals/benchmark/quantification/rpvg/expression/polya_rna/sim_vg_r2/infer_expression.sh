set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${GRAPH}_sim_vg_r2_${SIM}"
OUT_PREFIX="${QUANTER}_${MAPPER}_${TRANSCRIPTS}_sim_vg_r2_${SIM}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${GRAPH}/ . --recursive --exclude "*" --include "*.gam" --include "*.gamp" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt.ri . --no-progress

# Download transcript info
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.txt.gz" --exclude "*MT*" --exclude "*SCA*" --no-progress

# Concatenate transcript info
/usr/bin/time -v bash -c "zcat */*.txt.gz | head -n 1 > transcript_info.txt; zcat */*.txt.gz | grep -v ^Name >> transcript_info.txt; wc -l transcript_info.txt"

if [ "${MAPPER}" = "map" ] || [ "${MAPPER}" = "map_fast" ]; then 

	# Concatenate alignments
	/usr/bin/time -v bash -c "cat ${ALIGN_PREFIX}_*.gam > ${ALIGN_PREFIX}.gam"

	if [ "${QUANTER}" = "rpvg" ]; then 

		# Infer expression
		/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -u --score-not-qual -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"
	fi

elif [ "${MAPPER}" = "mpmap" ] || [ "${MAPPER}" = "mpmap_nosplice" ]; then
 
	# Concatenate alignments
	/usr/bin/time -v bash -c "cat ${ALIGN_PREFIX}_*.gamp > ${ALIGN_PREFIX}.gamp"

	if [ "${QUANTER}" = "rpvg" ]; then

		# Infer expression
		/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_trs" ]; then

		# Infer expression
		/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -i transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_gibbs" ]; then

		# Infer expression
		/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -n 1000 -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_gam" ]; then

		# Convert to gamp to gam
		/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"

		# Infer expression        
		/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -u -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"	
	fi
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"

if [ "${QUANTER}" != "rpvg_trs" ]; then

	# Compress expression values
	/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}_joint.txt"
fi

# Upload expression
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rpvg/expression/polya_rna/sim_vg_r2/${SIM}/${QUANTER}/${TRANSCRIPTS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
