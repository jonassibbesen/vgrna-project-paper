set -e

# Set transcripts (FASTA) prefix
TRANSCRIPTS_PREFIX="1kg_EURnonCEU_af002_gencode100"

# Set output name prefix
OUT_PREFIX="rsem_index_1kg_EURnonCEU_af002_gencode100"

# Set number of threads
CPU=1

# Construct rsem bowtie2 index
/usr/bin/time -v bash -c "rsem-prepare-reference -p ${CPU} --bowtie2 ${TRANSCRIPTS_PREFIX}.fa.gz ${OUT_PREFIX}"
