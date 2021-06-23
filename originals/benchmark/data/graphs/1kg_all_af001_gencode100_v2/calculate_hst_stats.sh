set -e

# Set file name prefixes
POPINFO_PREFIX="integrated_call_samples_v3.20130502.ALL"
GRAPHS_PREFIX="1kg_all_af001_gencode100_v2"
TRANSCRIPTS_PREFIX="${GRAPHS_PREFIX}_${CHR}"
OUT_PREFIX="${TRANSCRIPTS_PREFIX}_hst_stats"

# Download population info 
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/meta/1000g/${POPINFO_PREFIX}.panel . --no-progress

# Download haplotype-specific transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/${TRANSCRIPTS_PREFIX}.txt.gz . --no-progress

# Calculate haplotype-specific transcript statistics
/usr/bin/time -v bash -c "python3 /scripts/vgrna/haplotype_transcript_stats.py ${POPINFO_PREFIX}.panel <(zcat ${TRANSCRIPTS_PREFIX}.txt.gz) > ${OUT_PREFIX}.txt"

# Upload statistics
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
