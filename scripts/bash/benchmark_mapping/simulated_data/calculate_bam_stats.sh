set -e

# Set mapper name
MAPPER="mpmap"

# Set alignment (BAM) prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set simulated transcript alignments (BAM) prefix
TRANSCRIPTS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set simulated alignment positions (txt) prefix
SIM_PREFIX="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1"

# Set transcript annotation (gff) prefix
ANNOTATION_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set variants (vcf.gz) prefix
VARIANT_PREFIX="1kg_NA12878_exons"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Sort alignments by name
/usr/bin/time -v bash -c "samtools sort -n -O BAM ${ALIGN_PREFIX}.bam > ${ALIGN_PREFIX}_rsort.bam"

# Calculate alignment statistics
/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}_rsort.bam"

if [ "${MAPPER}" = "star" ] || [ "${MAPPER}" = "STAR" ]; then

	# Fix names for second reads in pair
	/usr/bin/time -v bash -c "samtools merge -O BAM ${ALIGN_PREFIX}_rsort_fixed.bam <(samtools view -h -f 128 ${ALIGN_PREFIX}_rsort.bam | awk -v OFS='\t' '{gsub("'"_1$"'", "'"_2"'", "'$1'")} 1' | samtools view -O BAM) <(samtools view -h -F 128 ${ALIGN_PREFIX}_rsort.bam | samtools view -O BAM)"

	ALIGN_PREFIX="${ALIGN_PREFIX}_rsort_fixed"

else

	ALIGN_PREFIX="${ALIGN_PREFIX}_rsort"
fi	

# Calculate overlap statistics
/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}.bam ${TRANSCRIPTS_PREFIX}.bam <(zcat ${SIM_PREFIX}.txt.gz) 3 ${ANNOTATION_PREFIX}.gff null ${VARIANT_PREFIX}.vcf.gz NA12878 null > ${OUT_PREFIX}_ovl3_vg.txt; gzip ${OUT_PREFIX}_ovl3_vg.txt"
