set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_real_${REAL}"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/real/${REAL}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download gencode exons
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full_exons.bed gencode_exons.bed --no-progress

# Calculate alignment statistics
/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}.bam"

# Calculate exon read overlap
/usr/bin/time -v bash -c "calc_read_regions_overlap_stats ${ALIGN_PREFIX}.bam gencode_exons.bed > ${OUT_PREFIX}_exon_ovl_gc.txt; gzip ${OUT_PREFIX}_exon_ovl_gc.txt"

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam gencode_exons.bed > ${OUT_PREFIX}_exon_cov_gc.txt; gzip ${OUT_PREFIX}_exon_cov_gc.txt"

if [ "${TRANSCRIPTS}" != "" ]; then

	# Download transcript alignments exons
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.bed" --no-progress

	# Calculate exon read overlap
	/usr/bin/time -v bash -c "calc_read_regions_overlap_stats ${ALIGN_PREFIX}.bam ${TRANSCRIPTS}_mq0.bed > ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq0.txt; gzip ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq0.txt"

	# Calculate exon read overlap
	/usr/bin/time -v bash -c "calc_read_regions_overlap_stats ${ALIGN_PREFIX}.bam ${TRANSCRIPTS}_mq30.bed > ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq30.txt; gzip ${OUT_PREFIX}_exon_ovl_${TRANSCRIPTS}_mq30.txt"


	# Calculate exon read coverage
	/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${TRANSCRIPTS}_mq0.bed > ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq0.txt; gzip ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq0.txt"

	# Calculate exon read coverage
	/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${TRANSCRIPTS}_mq30.bed > ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq30.txt; gzip ${OUT_PREFIX}_exon_cov_${TRANSCRIPTS}_mq30.txt"
fi

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/stats/polya_rna/real/${REAL}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
