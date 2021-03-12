set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_${SIM}"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

for i in $(seq 1 2); do

	# Filter secondary alignments
	/usr/bin/time -v bash -c "samtools view -F 256 -b ${OUT_PREFIX}_h${i}.bam > ${OUT_PREFIX}_h${i}_primary.bam"

	# Inject bam to gam
	/usr/bin/time -v bash -c "vg inject -t ${CPU} -x ${GRAPH}.xg ${OUT_PREFIX}_h${i}_primary.bam | vg view -a - | sed 's/\/1\",/\",/g' | sed 's/\/2\",/\",/g' | vg view -a -G -J - > ${OUT_PREFIX}_h${i}.gam"
done

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/sim_vg/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*.gam" --no-progress
