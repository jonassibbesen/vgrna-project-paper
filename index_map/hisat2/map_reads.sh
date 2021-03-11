set -e

# Set read files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_2.fq.gz"

# Set index prefix
INDEX_PREFIX="1kg_nonCEU_af001_gencode100_index"

# Set output name prefix
OUT_PREFIX="hisat2_1kg_nonCEU_af001_gencode100_sim_vg_SRR1153470"

# Set number of threads
CPU=1

# Map reads
/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t -x ${INDEX_PREFIX} -1 ${READ_1} -2 ${READ_2} -S ${OUT_PREFIX}.sam"

# Compress alignments
/usr/bin/time -v bash -c "samtools view -O BAM --threads ${CPU} ${OUT_PREFIX}.sam > ${OUT_PREFIX}.bam; rm ${OUT_PREFIX}.sam"

# Sort and index alignments
/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}.bam > ${OUT_PREFIX}_sort.bam; mv ${OUT_PREFIX}_sort.bam ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam"
