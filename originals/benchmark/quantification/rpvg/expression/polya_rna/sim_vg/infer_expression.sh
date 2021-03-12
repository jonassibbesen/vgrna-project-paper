set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_${SIM}"
OUT_PREFIX="${QUANTER}_${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.gam" --include "*.gamp" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/${REF}.xg . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/${REF}.gbwt . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/${REF}.gbwt.ri . --no-progress

# Download transcript info
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Concatenate transcript info
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep -v ^Name > transcript_info.txt; wc -l transcript_info.txt"

if [ "${MAPPER}" = "map" ]; then 

	# Concatenate alignments
	/usr/bin/time -v bash -c "cat ${ALIGN_PREFIX}_*.gam > ${ALIGN_PREFIX}.gam"

	if [ "${QUANTER}" = "rpvg" ]; then 

		# Infer expression
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${REF}.xg -p ${REF}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"
	fi

elif [ "${MAPPER}" = "mpmap" ] || [ "${MAPPER}" = "mpmap_nosplice" ]; then
 
	# Concatenate alignments
	/usr/bin/time -v bash -c "cat ${ALIGN_PREFIX}_*.gamp > ${ALIGN_PREFIX}.gamp"

	if [ "${QUANTER}" = "rpvg" ]; then

		# Infer expression
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -i haplotype-transcripts -g ${REF}.xg -p ${REF}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

	elif [ "${QUANTER}" = "rpvg_gam" ]; then

		# Convert to gamp to gam
		/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"

		# Infer expression        
		/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${REF}.xg -p ${REF}.gbwt -a ${ALIGN_PREFIX}.gam -f transcript_info.txt -o ${OUT_PREFIX}"	
	fi
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rpvg/expression/polya_rna/sim_vg/${SIM}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
