set -e

# Set population info prefixes
POPINFO_PREFIX="integrated_call_samples_v3.20130502.ALL"

# Set trancript info (txt) prefix
INFO_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_hst_stats"

# Set number of threads
CPU=1

# Calculate haplotype-specific transcript statistics
/usr/bin/time -v bash -c "python3 /scripts/vgrna/haplotype_transcript_stats.py ${POPINFO_PREFIX}.panel <(zcat ${INFO_PREFIX}.txt.gz) > ${OUT_PREFIX}.txt"
