set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"
OUT_PREFIX="rsem_index_${REF}"

# Download genome
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa.fai . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.fa.gz" --exclude "*MT*" --exclude "*SCA*" --no-progress

# Combine transcripts
/usr/bin/time -v bash -c "zcat */*fa.gz > ${REF}.fa; cat ${REF}.fa | grep '"'>'"' | wc -l"

# Construct rsem bowtie2 index
/usr/bin/time -v bash -c "rsem-prepare-reference -p ${CPU} --bowtie2 ${REF}.fa ${OUT_PREFIX}"

# Upload index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/quantification/rsem/indexes/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
