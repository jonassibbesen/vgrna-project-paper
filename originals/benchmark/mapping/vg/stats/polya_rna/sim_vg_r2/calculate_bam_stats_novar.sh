set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_r2_${SIM}"
SIM_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg_r1"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download transcript alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.bam . --no-progress

# Download simulated alignment positions
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg_r1/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878_exons/ . --recursive --exclude "*" --include "*.vcf.gz" --include "*.vcf.gz.tbi" --no-progress

if [ "${REF}" = "1kg_nonCEU_af001_gencode80" ]; then 

	# Download transcript annotation
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full_subset80.gff3 transcripts.gff3 --no-progress

else

	# Download transcript annotation
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full.gff3 transcripts.gff3 --no-progress
fi

# Combine variants
/usr/bin/time -v bash -c 'bcftools concat -O v -f <(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/1kg_NA12878_exons_${i}.vcf.gz; done) > 1kg_NA12878_exons.vcf; grep -v -P "^#" 1kg_NA12878_exons.vcf | wc -l; bgzip 1kg_NA12878_exons.vcf; tabix 1kg_NA12878_exons.vcf.gz'

for i in $(seq 1 2); do

	# Sort alignments by name
	/usr/bin/time -v bash -c "samtools sort -n -O BAM ${ALIGN_PREFIX}_h${i}.bam > ${ALIGN_PREFIX}_h${i}_rsort.bam"

	# Calculate alignment statistics
	/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}_h${i}_rsort.bam"

	# Calculate overlap statistics
	/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}_h${i}_rsort.bam ${TRANSCRIPTS}.bam <(zcat ${SIM_PREFIX}_h${i}.txt.gz) 3 transcripts.gff3 null null null > ${OUT_PREFIX}_bam_ovl3_h${i}.txt; gzip ${OUT_PREFIX}_bam_ovl3_h${i}.txt"
done

# Calculate allele coverage
/usr/bin/time -v bash -c "calc_allele_read_coverage ${ALIGN_PREFIX}_h1.bam ${ALIGN_PREFIX}_h2.bam <(zcat 1kg_NA12878_exons.vcf.gz) > ${OUT_PREFIX}_allele_cov.txt; gzip ${OUT_PREFIX}_allele_cov.txt"	 

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/stats/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
