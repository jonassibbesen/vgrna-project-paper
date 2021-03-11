set -e

# Set mapper name
MAPPER="mpmap"

# Set alignment prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set simulated transcript alignments prefix
TRANSCRIPTS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set simulated alignment positions prefix
SIM_PREFIX="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878_exons/ . --recursive --exclude "*" --include "*.vcf.gz" --include "*.vcf.gz.tbi" --no-progress

# Calculate alignment statistics
/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}.bam"

if [ "${MAPPER}" = "star" ]; then

	# Fix names for second reads in pair
	/usr/bin/time -v bash -c "samtools merge -O BAM ${ALIGN_PREFIX}_fixed.bam <(samtools view -h -f 128 ${ALIGN_PREFIX}.bam | awk -v OFS='\t' '{gsub("'"_1$"'", "'"_2"'", "'$1'")} 1' | samtools view -O BAM) <(samtools view -h -F 128 ${ALIGN_PREFIX}.bam | samtools view -O BAM)"

	ALIGN_PREFIX="${ALIGN_PREFIX}_fixed"
fi	

# Calculate overlap statistics
/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}.bam ${TRANSCRIPTS_PREFIX}.bam <(zcat ${SIM_PREFIX}.txt.gz) 3 > ${OUT_PREFIX}_ovl3_vg.txt; gzip ${OUT_PREFIX}_ovl3_vg.txt"
