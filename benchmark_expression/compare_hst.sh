set -e

# Set simulated transcript (FASTA) prefix
SIM_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set transcript (FASTA) prefix
TRANSCRIPTS_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set number of threads
CPU=1

# Compare hst sequences
/usr/bin/time -v bash -c "python3 /scripts/vgrna/compare_hst_sequences.py <(zcat ${SIM_PREFIX}.fa.gz) <(zcat ${TRANSCRIPTS_PREFIX}.fa.gz) ${OUT_PREFIX}.txt; wc -l ${OUT_PREFIX}.txt"
