set -e

# Set file name prefixes 
OUT_PREFIX="kallisto_index_${REF}"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.fa.gz" --no-progress

# Combine transcripts
/usr/bin/time -v bash -c "cat */*fa.gz > ${REF}.fa.gz; zcat ${REF}.fa.gz | grep '"'>'"' | wc -l"

# Construct kallisto index
/usr/bin/time -v bash -c "kallisto index -i ${OUT_PREFIX}.idx ${REF}.fa.gz"

# Upload index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/kallisto/indexes/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
