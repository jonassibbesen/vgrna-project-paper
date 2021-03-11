set -e

# Set expression profile prefix
PROFILE_PREFIX="1kg_NA12878_gencode100_ENCSR000AED_rep1_rsem"

# Set output name prefix
OUT_PREFIX="1kg_NA12878_gencode100_ENCSR000AED_rep1_uni"

# Set number of threads
CPU=1

# Average expression
/usr/bin/time -v bash -c "wc -l ${PROFILE_PREFIX}.isoforms.results; python3 /scripts/vgrna/uniform_diploid_rsem_expression.py ${PROFILE_PREFIX}.isoforms.results ${OUT_PREFIX}.isoforms.results; wc -l ${OUT_PREFIX}.isoforms.results"
