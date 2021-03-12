set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="1kg_NA12878_gencode100"
READS_PREFIX="ENCSR000AED_rep1"
EXPRESSION_PREFIX="${READS_PREFIX}_uni"
PROFILE_PREFIX="${TRANSCRIPTS_PREFIX}_${READS_PREFIX}_rsem"
OUT_PREFIX="${TRANSCRIPTS_PREFIX}_${EXPRESSION_PREFIX}"

# Download expression profile
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/${TRANSCRIPTS_PREFIX}/${READS_PREFIX}/${PROFILE_PREFIX}.isoforms.results . --no-progress

# Average expression
/usr/bin/time -v bash -c "wc -l ${PROFILE_PREFIX}.isoforms.results; python3 /scripts/vgrna/uniform_diploid_rsem_expression.py ${PROFILE_PREFIX}.isoforms.results ${OUT_PREFIX}.isoforms.results; wc -l ${OUT_PREFIX}.isoforms.results"

# Upload uniform expression
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/${TRANSCRIPTS_PREFIX}/${EXPRESSION_PREFIX}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
