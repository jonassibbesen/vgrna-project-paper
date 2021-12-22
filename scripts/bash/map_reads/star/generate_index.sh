set -e

# Set genome (FASTA) prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set transcripts (GTF) prefix
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set output name prefix
OUT_PREFIX="gencode100_index"

# Set number of threads
CPU=1

# Construct STAR index
/usr/bin/time -v bash -c "mkdir index; STAR --runThreadN ${CPU} --runMode genomeGenerate --genomeDir ./${OUT_PREFIX}/ --genomeFastaFiles ${GENOME_PREFIX}.fa --sjdbGTFfile ${TRANSCRIPTS_PREFIX}.gtf --sjdbOverhang 150"
