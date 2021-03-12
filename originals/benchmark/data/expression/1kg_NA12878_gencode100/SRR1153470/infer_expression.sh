set -e

# Set file name prefixes
TRANSCRIPTS_PREFIX="1kg_NA12878_gencode100"
READS_PREFIX="SRR1153470"
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
OUT_PREFIX="${TRANSCRIPTS_PREFIX}_${READS_PREFIX}_rsem"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/reads/illumina/SRP036136/${READS_PREFIX}_1.fastq.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/reads/illumina/SRP036136/${READS_PREFIX}_2.fastq.gz . --no-progress

# Download haplotype-specific transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.fa.gz" --include "*.txt.gz" --no-progress

# Extract haplotype-specific transcripts
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep NA12878 | cut -f1 | grep -v Name > haplotype_transcripts.txt; seqtk subseq <(zcat */*fa.gz) haplotype_transcripts.txt > ${TRANSCRIPTS_PREFIX}.fa; grep '"'>'"' ${TRANSCRIPTS_PREFIX}.fa | wc -l"

# Generate transcriptome reference and bowtie2 index
/usr/bin/time -v rsem-prepare-reference -p ${CPU} --bowtie2 ${TRANSCRIPTS_PREFIX}.fa ${TRANSCRIPTS_PREFIX} 

# Infer expression
/usr/bin/time -v rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --paired-end ${READS_PREFIX}_1.fastq.gz ${READS_PREFIX}_2.fastq.gz ${TRANSCRIPTS_PREFIX} ${OUT_PREFIX}

# Upload expression estimates and statistics
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/${TRANSCRIPTS_PREFIX}/${READS_PREFIX}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.bam"  --no-progress
