set -e

# Set file name prefixes
ALIGN_PREFIX="star_alleleseq_diploid_gencode100_real_r1"
OUT_PREFIX="${MAPPER}_${REF}_real_r2_${REAL}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jeizenga/alleleseq/ . --recursive --exclude "*" --include "${ALIGN_PREFIX}*${REAL}*${READS}*.bam" --no-progress

# Download gencode exons
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full_exons.bed gencode_exons.bed --no-progress

ALIGN_FILE=$(ls ${ALIGN_PREFIX}*)

# Sort and index alignments 
/usr/bin/time -v bash -c "samtools sort -o align_sort.bam ${ALIGN_FILE}; samtools index align_sort.bam"

# Calculate alignment statistics
/usr/bin/time -v bash -c "samtools flagstat align_sort.bam"

# Calculate exon read overlap
/usr/bin/time -v bash -c "calc_read_regions_overlap_stats align_sort.bam gencode_exons.bed > ${OUT_PREFIX}_exon_ovl_gc.txt; gzip ${OUT_PREFIX}_exon_ovl_gc.txt"

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage align_sort.bam gencode_exons.bed > ${OUT_PREFIX}_exon_cov_gc.txt; gzip ${OUT_PREFIX}_exon_cov_gc.txt"

if [ "${TRANSCRIPTS}" != "" ]; then

	# Download transcript alignments exons
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.bed" --no-progress

	# Calculate exon read overlap
	/usr/bin/time -v bash -c "calc_read_regions_overlap_stats align_sort.bam ${TRANSCRIPTS}_mq30.bed > ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq30.txt; gzip ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq30.txt"

	# Calculate exon read coverage
	/usr/bin/time -v bash -c "calc_exon_read_coverage align_sort.bam ${TRANSCRIPTS}_mq30.bed > ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq30.txt; gzip ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq30.txt"
fi

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/stats/polya_rna/real_r2/${REAL}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
