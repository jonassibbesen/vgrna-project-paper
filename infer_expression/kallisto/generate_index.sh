set -e

# Set transcripts prefix
TRANSCRIPTS_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="kallisto_index_1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

# Construct kallisto index
/usr/bin/time -v bash -c "kallisto index -i ${OUT_PREFIX}.idx ${TRANSCRIPTS_PREFIX}.fa.gz"
