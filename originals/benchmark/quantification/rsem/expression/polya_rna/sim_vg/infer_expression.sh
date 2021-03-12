set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg"
OUT_PREFIX="${QUANTER}_${REF}_sim_vg_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/indexes/${REF}/ . --recursive --no-progress

# Remove mate suffix 
/usr/bin/time -v bash -c "zcat ${READS_PREFIX}_h1_1.fq.gz ${READS_PREFIX}_h2_1.fq.gz | sed -e 's/_1$//g' | gzip > ${READS_PREFIX}_1.fq.gz; zcat ${READS_PREFIX}_h1_2.fq.gz ${READS_PREFIX}_h2_2.fq.gz | sed -e 's/_2$//g' | gzip > ${READS_PREFIX}_2.fq.gz"

# Use default RSEM
if [ "${QUANTER}" = "rsem" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --no-bam-output --paired-end ${READS_PREFIX}_1.fq.gz ${READS_PREFIX}_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"

# Use default RSEM with more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_k1k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 1000 --seed ${SEED} --no-bam-output --paired-end ${READS_PREFIX}_1.fq.gz ${READS_PREFIX}_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"

# Use default RSEM with many more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_k2k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 2000 --seed ${SEED} --no-bam-output --paired-end ${READS_PREFIX}_1.fq.gz ${READS_PREFIX}_2.fq.gz rsem_index_${REF} ${OUT_PREFIX}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.genes.results; gzip ${OUT_PREFIX}.isoforms.results"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/expression/polya_rna/sim_vg/${SIM}/${QUANTER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
