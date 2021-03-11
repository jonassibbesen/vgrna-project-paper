set -e

# Set input name prefixes
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
VARIANTS_PREFIX="1kg_all_exons_${CHR}"
GRAPH_PREFIX="1kg_all_af001_gencode100_${CHR}"

# Set output name prefixes
OUT_PREFIX="1kg_all_af001_gencode100_${CHR}"

# Set number of threads
CPU=1

# Create gbwt index of all haplotypes
/usr/bin/time -v bash -c "vg index -p -t ${CPU} -G ${VARIANTS_PREFIX}.gbwt -v ${VARIANTS_PREFIX}.vcf.gz ${GRAPH_PREFIX}.pg"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${GRAPH_PREFIX}.pg"

# Create haplotype-specific transcripts and update graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -r -g -n ${CHR}.gtf -l ${VARIANTS_PREFIX}.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg"

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"
