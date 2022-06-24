set -e

# Set file name prefixes 
OUT_PREFIX="rsem_index_${REF}_mt_sca"

# Download main transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.fa.gz" --no-progress

# Download MT and SCA transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}_mt_sca/ . --recursive --exclude "*" --include "*.fa.gz" --no-progress

# Combine transcripts
/usr/bin/time -v bash -c "zcat */*fa.gz > ${REF}.fa; cat ${REF}.fa | grep '"'>'"' | wc -l"

# Construct rsem bowtie2 index
/usr/bin/time -v bash -c "rsem-prepare-reference -p ${CPU} --bowtie2 ${REF}.fa ${OUT_PREFIX}"

# Upload index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/indexes/${REF}_mt_sca/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
