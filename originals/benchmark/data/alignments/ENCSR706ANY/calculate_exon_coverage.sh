set -e

# Set file name prefixes
ENCODE_ID="ENCSR706ANY"
ALIGN_PREFIX="${ENCODE_ID}_mq${MAPQ}"
OUT_PREFIX="${ALIGN_PREFIX}_exon_cov"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${ENCODE_ID}/${ALIGN_PREFIX}.bam . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${ENCODE_ID}/${ALIGN_PREFIX}.bam.bai . --no-progress

# Download alignment exons
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${ENCODE_ID}/${ALIGN_PREFIX}.bed . --no-progress

# Download gencode exons
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full_exons.bed exons_gc.bed --no-progress

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam ${ALIGN_PREFIX}.bed > ${OUT_PREFIX}_bam.txt; gzip ${OUT_PREFIX}_bam.txt"

# Calculate exon read coverage
/usr/bin/time -v bash -c "calc_exon_read_coverage ${ALIGN_PREFIX}.bam exons_gc.bed > ${OUT_PREFIX}_gc.txt; gzip ${OUT_PREFIX}_gc.txt"

# Upload exon coverage 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/alignments/${ENCODE_ID}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
