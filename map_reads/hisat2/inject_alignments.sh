set -e

# Set alignment (BAM) prefix
ALIGN_PREFIX="hisat2_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set graph (XG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="hisat2_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Filter secondary alignments
/usr/bin/time -v bash -c "samtools view -F 256 -b ${OUT_PREFIX}.bam > ${OUT_PREFIX}_primary.bam"

# Inject bam to gam
/usr/bin/time -v bash -c "vg inject -t ${CPU} -x ${GRAPH_PREFIX}.xg ${OUT_PREFIX}_primary.bam | vg view -a - | sed 's/\/1\",/\",/g' | sed 's/\/2\",/\",/g' | vg view -a -G -J - > ${OUT_PREFIX}.gam"
