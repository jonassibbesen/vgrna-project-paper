set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_${SIM}"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.gam" --include "*.gamp" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/${REF}.xg . --no-progress

for i in $(seq 1 2); do

	# Extract reference path list
	/usr/bin/time -v bash -c "vg paths -L -x ${REF}.xg | grep -v ENST > reference_paths.txt; wc -l reference_paths.txt"

	if [ -f "${ALIGN_PREFIX}_h${i}.gamp" ]; then

		# Surject gam to bam
		/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -A -b -m -F reference_paths.txt -x ${REF}.xg ${ALIGN_PREFIX}_h${i}.gamp > ${OUT_PREFIX}_h${i}.bam"

	else 

		# Surject gam to bam
		/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -A -b -F reference_paths.txt -x ${REF}.xg ${ALIGN_PREFIX}_h${i}.gam > ${OUT_PREFIX}_h${i}.bam"
	fi

	# Sort, compress and index alignments
	/usr/bin/time -v bash -c "samtools sort -O BAM --threads ${CPU} ${OUT_PREFIX}_h${i}.bam > ${OUT_PREFIX}_h${i}_sort.bam; mv ${OUT_PREFIX}_h${i}_sort.bam ${OUT_PREFIX}_h${i}.bam; samtools index ${OUT_PREFIX}_h${i}.bam"
done 

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*.bam*" --no-progress
