set -e

# Set inference mode. 
# 	"kallisto": Use default
#	"kallisto_strand": Use strand-specific
#	"kallisto_strand_bias": Use strand-specific with bias-correction
MODE="kallisto"

# Set rng seed
SEED=622797

# Set read files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_2.fq.gz"

# Set index prefix
INDEX_PREFIX="kallisto_index_1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="kallisto_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

# Use Kallisto
if [ "${QUANTER}" = "kallisto" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} -i ${INDEX_PREFIX}.idx -o ${OUT_PREFIX} ${READ_1} ${READ_2}"

# Use stranded Kallisto
elif [ "${QUANTER}" = "kallisto_strand" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} --rf-stranded -i ${INDEX_PREFIX}.idx -o ${OUT_PREFIX} ${READ_1} ${READ_2}"

# Use Salmon with bias correction
elif [ "${QUANTER}" = "kallisto_strand_bias" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "kallisto quant -t ${CPU} --seed ${SEED} --rf-stranded --bias -i ${INDEX_PREFIX}.idx -o ${OUT_PREFIX} ${READ_1} ${READ_2}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/abundance.tsv"
