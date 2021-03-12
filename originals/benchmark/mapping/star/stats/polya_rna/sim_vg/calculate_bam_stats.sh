set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_${SIM}"
SIM_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download transcript alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.bam . --no-progress

# Download simulated alignment positions
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878_exons/ . --recursive --exclude "*" --include "*.vcf.gz" --include "*.vcf.gz.tbi" --no-progress

for i in $(seq 1 2); do

	# Calculate alignment statistics
	/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}_h${i}.bam"
	
	# Fix names for second reads in pair
	/usr/bin/time -v bash -c "samtools merge -O BAM ${ALIGN_PREFIX}_h${i}_fixed.bam <(samtools view -h -f 128 ${ALIGN_PREFIX}_h${i}.bam | awk -v OFS='\t' '{gsub("'"_1$"'", "'"_2"'", "'$1'")} 1' | samtools view -O BAM) <(samtools view -h -F 128 ${ALIGN_PREFIX}_h${i}.bam | samtools view -O BAM)"

	# Calculate overlap statistics
	/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}_h${i}_fixed.bam ${TRANSCRIPTS}.bam <(zcat ${SIM_PREFIX}_h${i}.txt.gz) 0 > ${OUT_PREFIX}_ovl0_vg_h${i}.txt; gzip ${OUT_PREFIX}_ovl0_vg_h${i}.txt"

	# Calculate overlap statistics
	/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}_h${i}_fixed.bam ${TRANSCRIPTS}.bam <(zcat ${SIM_PREFIX}_h${i}.txt.gz) 3 > ${OUT_PREFIX}_ovl3_vg_h${i}.txt; gzip ${OUT_PREFIX}_ovl3_vg_h${i}.txt"

done

# Combine variants
/usr/bin/time -v bash -c 'bcftools concat -O v -f <(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/1kg_NA12878_exons_${i}.vcf.gz; done) > 1kg_NA12878_exons.vcf'

# Calculate allele coverage
/usr/bin/time -v bash -c "calc_allele_read_coverage ${ALIGN_PREFIX}_h1.bam ${ALIGN_PREFIX}_h2.bam 1kg_NA12878_exons.vcf > ${OUT_PREFIX}_allele_cov.txt; gzip ${OUT_PREFIX}_allele_cov.txt"	 

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/stats/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
