set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg"
OUT_PREFIX="${QUANTER}_${REF}_sim_vg_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/salmon/indexes/${REF}/ . --recursive --no-progress

# Concatenate reads
/usr/bin/time -v bash -c "cat ${READS_PREFIX}_h1_1.fq.gz ${READS_PREFIX}_h2_1.fq.gz > ${READS_PREFIX}_1.fq.gz; cat ${READS_PREFIX}_h1_2.fq.gz ${READS_PREFIX}_h2_2.fq.gz > ${READS_PREFIX}_2.fq.gz"

# Use default Salmon
if [ "${QUANTER}" = "salmon" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with additional read mappings
elif [ "${QUANTER}" = "salmon_w5k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A -w 5000 -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with gibbs sampling (100 samples)
elif [ "${QUANTER}" = "salmon_gibbs100" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --numGibbsSamples 100 -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with gibbs sampling (1000 samples)
elif [ "${QUANTER}" = "salmon_gibbs1000" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --numGibbsSamples 1000 -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with traditional EM algorithm
elif [ "${QUANTER}" = "salmon_em" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --useEM -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with a transcript-level VB prior of 0.1
elif [ "${QUANTER}" = "salmon_vbp01" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --perTranscriptPrior --vbPrior 0.1 -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"

# Use Salmon with a transcript-level VB prior of 1
elif [ "${QUANTER}" = "salmon_vbp1" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --perTranscriptPrior --vbPrior 1 -i salmon_index_${REF} -o ${OUT_PREFIX} -1 ${READS_PREFIX}_1.fq.gz -2 ${READS_PREFIX}_2.fq.gz"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/quant.sf"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/salmon/expression/polya_rna/sim_vg/${SIM}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
