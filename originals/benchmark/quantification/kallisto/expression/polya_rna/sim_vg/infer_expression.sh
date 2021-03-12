set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg"
OUT_PREFIX="${QUANTER}_${REF}_sim_vg_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/kallisto/indexes/${REF}/ . --recursive --no-progress

# Concatenate reads
/usr/bin/time -v bash -c "cat ${READS_PREFIX}_h1_1.fq.gz ${READS_PREFIX}_h2_1.fq.gz > ${READS_PREFIX}_1.fq.gz; cat ${READS_PREFIX}_h1_2.fq.gz ${READS_PREFIX}_h2_2.fq.gz > ${READS_PREFIX}_2.fq.gz"

# Use default Kallisto
if [ "${QUANTER}" = "kallisto" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} -i kallisto_index_${REF}.idx -o ${OUT_PREFIX} ${READS_PREFIX}_1.fq.gz ${READS_PREFIX}_2.fq.gz"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/abundance.tsv"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/kallisto/expression/polya_rna/sim_vg/${SIM}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
