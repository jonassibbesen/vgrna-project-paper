set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_r2_${SIM}"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.bam" --include "*.bam.bai" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

for i in $(seq 1 2); do

	# Inject bam to gam
	/usr/bin/time -v bash -c "vg inject -t ${CPU} -x ${GRAPH}.xg ${OUT_PREFIX}_h${i}.bam > align_tmp_h${i}.gam"

	# Fix read names
	/usr/bin/time -v bash -c "vg view -a align_tmp_h${i}.gam | sed 's/\/1\",/\",/g' | sed 's/_1\/2\",/_2\",/g' | vg view -a -G -J - > ${OUT_PREFIX}_h${i}.gam"
done 

# Upload statistics
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*.gam" --no-progress
