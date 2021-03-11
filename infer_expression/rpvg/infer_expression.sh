set -e

# Set inference mode. 
# 	"rpvg": Use default
#	"rpvg_gam": Use single-path alignments
#	"rpvg_strand": Use strand-specific
#	"rpvg_strand_gam": Use strand-specific and single-path alignments
QUANTER="rpvg"

# Set rng seed
SEED=622797

# Set alignment prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_genes_sim_vg_ENCSR000AED_rep1_uni"

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100_genes"

# Set trancript index prefix
INDEX_PREFIX="1kg_nonCEU_af001_gencode100_genes"

# Set trancript info prefix
INFO_PREFIX="1kg_nonCEU_af001_gencode100_genes"

# Set output name prefix
OUT_PREFIX="rpvg_1kg_nonCEU_af001_gencode100_genes_sim_vg_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

if [ "${QUANTER}" = "rpvg" ]; then 

	# Infer expression
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gamp -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"

elif [ "${QUANTER}" = "rpvg_gam" ]; then

	# Convert to gamp to gam
	/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"

	# Infer expression        
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gam -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"	

elif [ "${QUANTER}" = "rpvg_strand" ]; then 

	# Infer expression
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gamp -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"

elif [ "${QUANTER}" = "rpvg_strand_gam" ]; then

	# Convert to gamp to gam
	/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"

	# Infer expression        
	/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -u -m ${MEAN} -d ${SD} -i haplotype-transcripts -g ${GRAPH_PREFIX}.xg -p ${INDEX_PREFIX}.gbwt -a ${ALIGN_PREFIX}.gam -f <(zcat ${INFO_PREFIX}.txt.gz) -o ${OUT_PREFIX}"
fi

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"
