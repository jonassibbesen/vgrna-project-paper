set -e

# Set read (FASTQ) files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_2.fq.gz"

# Set index directory
INDEX_DIR="gencode100_index"

# Set output name prefix
OUT_PREFIX="star_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Map reads
/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir ${INDEX_DIR} --readFilesCommand zcat --readFilesIn ${READ_1} ${READ_2} --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_; mv ${OUT_PREFIX}_Aligned.out.bam ${OUT_PREFIX}.bam"

# Sort and index alignments
/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}.bam > ${OUT_PREFIX}_sort.bam; mv ${OUT_PREFIX}_sort.bam ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam"
