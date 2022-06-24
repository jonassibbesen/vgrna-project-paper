set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg_r1"
OUT_PREFIX="${MAPPER}_${REF}_sim_vg_r2_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg_r1/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/indexes/${REF}/ . --recursive --no-progress

for i in $(seq 1 2); do

	# Use default HISAT2
	if [ "${MAPPER}" = "hisat2" ]; then

		# Map reads
		/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t -x ${REF}_index -1 ${READS_PREFIX}_h${i}_1.fq.gz -2 ${READS_PREFIX}_h${i}_2.fq.gz -S ${OUT_PREFIX}_h${i}.sam"

	# Use HISAT2 with multimapping
	elif [ "${MAPPER}" = "hisat2_multi10" ]; then

		# Map reads
		/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t -k 10 --secondary -x ${REF}_index -1 ${READS_PREFIX}_h${i}_1.fq.gz -2 ${READS_PREFIX}_h${i}_2.fq.gz -S ${OUT_PREFIX}_h${i}.sam"
	fi

	# Compress alignments
	/usr/bin/time -v bash -c "samtools view -O BAM --threads ${CPU} ${OUT_PREFIX}_h${i}.sam > ${OUT_PREFIX}_h${i}.bam"

	# Sort and index alignments
	/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}_h${i}.bam > ${OUT_PREFIX}_h${i}_sort.bam; mv ${OUT_PREFIX}_h${i}_sort.bam ${OUT_PREFIX}_h${i}.bam; samtools index ${OUT_PREFIX}_h${i}.bam"
done

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.sam" --no-progress
