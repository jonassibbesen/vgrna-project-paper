set -e

# Set rng seed
SEED=828689658

# Set read (FASTQ) files
READ_1="ENCFF001REK.fastq.gz"
READ_2="ENCFF001REJ.fastq.gz"

# Set transcripts (FASTA) prefix
TRANSCRIPTS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set output name prefix
OUT_PREFIX="1kg_NA12878_gencode100_ENCSR000AED_rep1_rsem"

# Set number of threads
CPU=1

# Extract haplotype-specific transcripts
/usr/bin/time -v bash -c "zcat ${TRANSCRIPTS_PREFIX}.txt.gz | grep NA12878 | cut -f1 | grep -v Name > haplotype_transcripts.txt; seqtk subseq <(${TRANSCRIPTS_PREFIX}.fa.gz) haplotype_transcripts.txt > ${OUT_PREFIX}.fa; grep '"'>'"' ${OUT_PREFIX}.fa | wc -l; rm haplotype_transcripts.txt"

# Generate transcriptome reference and bowtie2 index
/usr/bin/time -v rsem-prepare-reference -p ${CPU} --bowtie2 ${OUT_PREFIX}.fa ${OUT_PREFIX} 

# Infer expression
/usr/bin/time -v rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --paired-end ${READ_1} ${READ_2} ${OUT_PREFIX} ${OUT_PREFIX}
