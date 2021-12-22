set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_rsem_r1"
OUT_PREFIX="${MAPPER}_${REF}_sim_rsem_r1_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/rsem_r1/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/indexes/${REF}/ . --recursive --no-progress

for i in $(seq 1 2); do

	# Use default STAR
	if [ "${MAPPER}" = "star" ]; then

		# Map reads
		/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir . --readFilesCommand zcat --readFilesIn ${READS_PREFIX}_h${i}_1.fq.gz ${READS_PREFIX}_h${i}_2.fq.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_h${i}_; mv ${OUT_PREFIX}_h${i}_Aligned.out.bam ${OUT_PREFIX}_h${i}.bam"

		# Sort and index alignments
		/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}_h${i}.bam > ${OUT_PREFIX}_h${i}_sort.bam; mv ${OUT_PREFIX}_h${i}_sort.bam ${OUT_PREFIX}_h${i}.bam; samtools index ${OUT_PREFIX}_h${i}.bam"
	fi
done

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/sim_rsem_r1/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.sam" --no-progress
