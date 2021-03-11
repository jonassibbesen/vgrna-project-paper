set -e

# Set alignment prefix
ALIGN_PREFIX="ENCSR706ANY_mq30"

# Set exons prefix
EXONS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full_exons"

# Set output name prefix
OUT_PREFIX="ENCSR706ANY_mq30_exon_cov"

# Set number of threads
CPU=1

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${ALIGN_PREFIX}.bed > ${OUT_PREFIX}_bam.txt; gzip ${OUT_PREFIX}_bam.txt"

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${EXONS_PREFIX}.bed > ${OUT_PREFIX}_gc.txt; gzip ${OUT_PREFIX}_gc.txt"
