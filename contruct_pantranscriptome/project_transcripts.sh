set -e

# Set projection mode. 
# 	"def": Default
# 	"redu": Keep redundant haplotype-specific transcripts
# 	"bidi": Create bi-directional pantranscriptome
MODE="def"

# Set transcripts (GTF) prefix
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set exon variants (VCF) prefix
VARIANTS_PREFIX="1kg_nonCEU_exons_${CHR}"

# Set graph PG prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100_${CHR}"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_${CHR}"

# Set number of threads
CPU=1

# Create gbwt index of all haplotypes
/usr/bin/time -v bash -c "vg index -p -t ${CPU} -G haplotypes.gbwt -v ${VARIANTS_PREFIX}.vcf.gz ${GRAPH_PREFIX}.pg"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${GRAPH_PREFIX}.pg"

if [ "${MODE}" = "def" ]; then 

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -r -n ${CHR}.gtf -l haplotypes.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg; rm haplotypes.gbwt"

elif [ "${MODE}" = "redu" ]; then

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -c -r -a -u -g -n ${CHR}.gtf -l haplotypes.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${OUT_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg; rm haplotypes.gbwt"

elif [ "${MODE}" = "bidi" ]; then

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -r -g -n ${CHR}.gtf -l haplotypes.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg; rm haplotypes.gbwt"
fi

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"
