set -e

# Set input name prefixes
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation"

# Set output name prefixes
OUT_PREFIX="${TRANSCRIPTS_PREFIX}_renamed_full"

# Set number of threads
CPU=1

# Filter for full-length transcripts and rename contigs
zcat ${TRANSCRIPTS_PREFIX}.gtf.gz | sed -e 's/^chrM/chrMT/g' | grep -v "mRNA_start_NF" | grep -v "mRNA_end_NF" | sed -e 's/^chr//g' > ${OUT_PREFIX}.gtf

# Select exons
grep -v "^#" ${OUT_PREFIX}.gtf | grep -P "\texon\t" | cut -f1,4,5 | awk -v OFS="\t" '{$2 -= 1}{print}' > ${OUT_PREFIX}_exons.txt

# Sort and merge exons
bedtools sort -i ${OUT_PREFIX}_exons.txt > ${OUT_PREFIX}_exons_sort.txt; bedtools merge -i ${OUT_PREFIX}_exons_sort.txt > ${OUT_PREFIX}_exons.bed; rm ${OUT_PREFIX}_exons*.txt
