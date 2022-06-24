set -e

# Set alignment (BAM) prefix
ALIGN_PREFIX="hisat2_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set graph (XG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="hisat2_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Inject bam to gam
/usr/bin/time -v bash -c "vg inject -t ${CPU} -x ${GRAPH_PREFIX}.xg ${OUT_PREFIX}.bam > ${OUT_PREFIX}.gam"

# Fix read names
/usr/bin/time -v bash -c "vg view -a ${OUT_PREFIX}.gam | sed 's/\/1\",/\",/g' | sed 's/\/2\",/\",/g' | vg view -a -G -J - > ${OUT_PREFIX}_fixed.gam"
