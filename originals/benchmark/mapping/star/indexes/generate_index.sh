set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Download genome
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa.fai . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS} . --no-progress

# Construct STAR index
/usr/bin/time -v bash -c "mkdir index; STAR --runThreadN ${CPU} --runMode genomeGenerate --genomeDir ./index/ --genomeFastaFiles ${GENOME_PREFIX}.fa --sjdbGTFfile ${TRANSCRIPTS} --sjdbOverhang 150"

# Upload STAR index
aws s3 sync ./index/ s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/indexes/${REF}/ --no-progress
