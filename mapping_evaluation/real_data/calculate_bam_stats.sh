set -e

# Set alignment prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_real_ENCSR000AED_rep1"

# Set exons prefix
EXONS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full_exons"

# Set regions prefix
REGIONS_PREFIX="ENCSR706ANY_mq30_exon_cov"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_real_ENCSR000AED_rep1"

# Set number of threads
CPU=1

# Calculate alignment statistics
/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}.bam"

# Calculate exon read overlap
/usr/bin/time -v bash -c "calc_read_regions_overlap_stats ${ALIGN_PREFIX}.bam ${EXONS_PREFIX}.bed > ${OUT_PREFIX}_exon_ovl_gc.txt; gzip ${OUT_PREFIX}_exon_ovl_gc.txt"

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${EXONS_PREFIX}.bed > ${OUT_PREFIX}_exon_cov_gc.txt; gzip ${OUT_PREFIX}_exon_cov_gc.txt"

if [ "${REGIONS_PREFIX}" != "" ]; then

	# Calculate exon read overlap
	/usr/bin/time -v bash -c "calc_read_regions_overlap_stats ${ALIGN_PREFIX}.bam ${REGIONS_PREFIX}.bed > ${OUT_PREFIX}_exon_ovl_${REGIONS_PREFIX}.txt; gzip ${OUT_PREFIX}_exon_ovl_${REGIONS_PREFIX}.txt"

	# Calculate exon read coverage
	/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${REGIONS_PREFIX}.bed > ${OUT_PREFIX}_exon_cov_${REGIONS_PREFIX}.txt; gzip ${OUT_PREFIX}_exon_cov_${REGIONS_PREFIX}.txt"
fi
