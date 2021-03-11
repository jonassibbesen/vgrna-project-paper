set -e

# Set alignment prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Extract reference path list
/usr/bin/time -v bash -c "vg paths -L -x ${GRAPH_PREFIX}.xg | grep -v ENST > reference_paths.txt; wc -l reference_paths.txt"

if [ -f "${ALIGN_PREFIX}.gamp" ]; then

	# Surject gamp to bam
	/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -A -b -m -F reference_paths.txt -x ${GRAPH_PREFIX}.xg ${ALIGN_PREFIX}.gamp > ${OUT_PREFIX}.bam"

else

	# Surject gam to bam
	/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -A -b -F reference_paths.txt -x ${GRAPH_PREFIX}.xg ${ALIGN_PREFIX}.gam > ${OUT_PREFIX}.bam"
fi

# Sort, compress and index alignments
/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}.bam > ${OUT_PREFIX}_sort.bam; mv ${OUT_PREFIX}_sort.bam ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam"
