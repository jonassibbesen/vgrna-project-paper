set -e

# Set inference mode. 
# 	"rsem": Use default
# 	"rsem_k1k": Use more Bowtie2 mappings
#	"rsem_strand_k1k": Use strand-specific and more Bowtie2 mappings
QUANTER="rsem"

# Set rng seed
SEED=622797

# Set read files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_2.fq.gz"

# Set index prefix
INDEX_PREFIX="rsem_index_1kg_EURnonCEU_af002_gencode100"

# Set output name prefix
OUT_PREFIX="rsem_1kg_EURnonCEU_af002_gencode100_sim_vg_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

# Use stranded RSEM
if [ "${QUANTER}" = "rsem_strand" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --no-bam-output --strandedness reverse --paired-end ${READ_1} ${READ_2} ${INDEX_PREFIX} ${OUT_PREFIX}"

# Use default RSEM with more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_k1k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 1000 --seed ${SEED} --no-bam-output --paired-end ${READ_1} ${READ_2} ${INDEX_PREFIX} ${OUT_PREFIX}"

# Use stranded RSEM with many more Bowtie2 mappings
elif [ "${QUANTER}" = "rsem_strand_k1k" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "rsem-calculate-expression -p ${CPU} --bowtie2 --bowtie2-k 1000 --seed ${SEED} --no-bam-output --strandedness reverse --paired-end ${READ_1} ${READ_2} ${INDEX_PREFIX} ${OUT_PREFIX}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.genes.results; gzip ${OUT_PREFIX}.isoforms.results"
