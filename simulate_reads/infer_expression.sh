set -e

# Set graph prefix
GRAPH_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set transcripts prefix
TRANSCRIPTS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set output name prefix
OUT_PREFIX="1kg_NA12878_gencode100_ENCSR000AED_rep1_rsem"

# Set number of threads
CPU=1

# Extract haplotype-specific transcripts
/usr/bin/time -v bash -c "zcat ${TRANSCRIPTS_PREFIX}.txt.gz | grep NA12878 | cut -f1 | grep -v Name > haplotype_transcripts.txt; seqtk subseq <(zcat */*fa.gz) haplotype_transcripts.txt > ${OUT_PREFIX}.fa; grep '"'>'"' ${OUT_PREFIX}.fa | wc -l; rm haplotype_transcripts.txt"

# Generate transcriptome reference and bowtie2 index
/usr/bin/time -v rsem-prepare-reference -p ${CPU} --bowtie2 ${OUT_PREFIX}.fa ${OUT_PREFIX} 

# Infer expression
/usr/bin/time -v rsem-calculate-expression -p ${CPU} --bowtie2 --seed ${SEED} --paired-end ENCFF001REK.fastq.gz ENCFF001REJ.fastq.gz ${OUT_PREFIX} ${OUT_PREFIX}
