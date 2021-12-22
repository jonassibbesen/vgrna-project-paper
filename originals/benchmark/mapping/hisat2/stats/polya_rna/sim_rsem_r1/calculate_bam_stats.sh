set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_rsem_r1_${SIM}"
SIM_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_rsem_r1"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/sim_rsem_r1/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download transcript alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.bam . --no-progress

# Download simulated transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/rsem_r1/ . --recursive --exclude "*" --include "*.isoforms.results" --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878_exons/ . --recursive --exclude "*" --include "*.vcf.gz" --include "*.vcf.gz.tbi" --no-progress

# Combine variants
/usr/bin/time -v bash -c 'bcftools concat -O v -f <(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/1kg_NA12878_exons_${i}.vcf.gz; done) > 1kg_NA12878_exons.vcf; grep -v -P "^#" 1kg_NA12878_exons.vcf | wc -l; bgzip 1kg_NA12878_exons.vcf; tabix 1kg_NA12878_exons.vcf.gz'

for i in $(seq 1 2); do

	# Calculate alignment statistics
	/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}_h${i}.bam"

	# Create vg sim position table
	/usr/bin/time -v bash -c "convert_rsem_sim_info ${ALIGN_PREFIX}_h${i}.bam ${SIM_PREFIX}_h${i}.sim.isoforms.results > ${ALIGN_PREFIX}_h${i}_pos.txt"

	# Calculate overlap statistics
	/usr/bin/time -v bash -c "calc_vg_benchmark_stats ${ALIGN_PREFIX}_h${i}.bam ${TRANSCRIPTS}.bam ${ALIGN_PREFIX}_h${i}_pos.txt 3 1kg_NA12878_exons.vcf.gz NA12878 > ${OUT_PREFIX}_ovl3_vg_h${i}.txt; gzip ${OUT_PREFIX}_ovl3_vg_h${i}.txt"
done

# Calculate allele coverage
/usr/bin/time -v bash -c "calc_allele_read_coverage ${ALIGN_PREFIX}_h1.bam ${ALIGN_PREFIX}_h2.bam <(zcat 1kg_NA12878_exons.vcf.gz) > ${OUT_PREFIX}_allele_cov.txt; gzip ${OUT_PREFIX}_allele_cov.txt"	 

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/stats/polya_rna/sim_rsem_r1/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
