set -e

# Set inference mode. 
# 	"rpvg": Use default
#	"rpvg_strand": Use strand-specific
MODE="rpvg"

# Set rng seed
SEED=622797

# Set alignment (GAM or GAMP) prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni"

# Set graph (XG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set trancript index (GBWT) prefix
INDEX_PREFIX="1kg_nonCEU_af001_gencode100"

# Set trancript info (txt) prefix
INFO_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="rpvg_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

if [ "${MODE}" = "rpvg" ]; then 

	# Infer expression
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gamp -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"

elif [ "${MODE}" = "rpvg_strand" ]; then 

	# Infer expression
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gamp -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"
