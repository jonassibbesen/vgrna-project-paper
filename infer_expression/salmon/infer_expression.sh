set -e

# Set inference mode. 
# 	"salmon": Use default
#	"salmon_bias": Use with bias-correction
MODE="salmon"

# Set read files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_2.fq.gz"

# Set index prefix
INDEX_PREFIX="salmon_index_1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="salmon_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

# Use default Salmon
if [ "${MODE}" = "salmon" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A -i ${INDEX_PREFIX} -o ${OUT_PREFIX} -1 ${READ_1} -2 ${READ_2}"

# Use Salmon with bias correction
elif [ "${MODE}" = "salmon_bias" ]; then

	# Infer expression
	/usr/bin/time -v bash -c "salmon quant -p ${CPU} -l A --seqBias --gcBias -i ${INDEX_PREFIX} -o ${OUT_PREFIX} -1 ${READ_1} -2 ${READ_2}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}/quant.sf"
