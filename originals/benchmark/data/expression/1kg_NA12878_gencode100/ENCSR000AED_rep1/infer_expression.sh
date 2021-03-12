set -e

# Set file name prefixes
TRANSCRIPTS_PREFIX="1kg_NA12878_gencode100"
READS_PREFIX="ENCSR000AED_rep1"
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
OUT_PREFIX="${TRANSCRIPTS_PREFIX}_${READS_PREFIX}_rsem"

# Download reads
aws s3 cp s3://encode-public/2013/06/13/c653a32e-e618-42b1-b8b8-b3b838847b97/ENCFF001REK.fastq.gz . --no-sign-request --no-progress
aws s3 cp s3://encode-public/2013/06/13/efa1a02d-6b43-4635-9ef8-d2d78c527839/ENCFF001REJ.fastq.gz . --no-sign-request --no-progress

# Download haplotype-specific transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.fa.gz" --include "*.txt.gz" --no-progress

# Extract haplotype-specific transcripts
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep NA12878 | cut -f1 | grep -v Name > haplotype_transcripts.txt; seqtk subseq <(zcat */*fa.gz) haplotype_transcripts.txt > ${TRANSCRIPTS_PREFIX}.fa; grep '"'>'"' ${TRANSCRIPTS_PREFIX}.fa | wc -l"

# Generate transcriptome reference and bowtie2 index
/usr/bin/time -v rsem-prepare-reference -p ${CPU} --bowtie2 ${TRANSCRIPTS_PREFIX}.fa ${TRANSCRIPTS_PREFIX} 

# Infer expression
/usr/bin/time -v rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --paired-end ENCFF001REK.fastq.gz ENCFF001REJ.fastq.gz ${TRANSCRIPTS_PREFIX} ${OUT_PREFIX}

# Upload expression estimates and statistics
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/${TRANSCRIPTS_PREFIX}/${READS_PREFIX}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.bam"  --no-progress
