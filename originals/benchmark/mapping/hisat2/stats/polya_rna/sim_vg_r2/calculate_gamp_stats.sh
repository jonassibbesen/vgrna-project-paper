set -e

# Set file name prefixes
ALIGN_PREFIX="${MAPPER}_${REF}_sim_vg_r2_${SIM}"
SIM_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_vg_r1"
OUT_PREFIX="${ALIGN_PREFIX}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/alignments/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ . --recursive --exclude "*" --include "*.gam" --include "*.gamp" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

# Download simulated alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/vg_r1/ . --recursive --exclude "*" --include "*.gam" --no-progress

for i in $(seq 1 2); do

	# Calculate alignment statistics
	/usr/bin/time -v bash -c "vg stats -a ${ALIGN_PREFIX}_h${i}.gam"

	# Calculate distance statistics
	/usr/bin/time -v bash -c "vg gampcompare -t ${CPU} -d -a ${MAPPER} -G ${GRAPH}.xg ${ALIGN_PREFIX}_h${i}.gam ${SIM_PREFIX}_h${i}.gam | grep -v ^distance | sort -nr -k4 | sort -ns -k1 | sort -su -k6 | cut -f1-5 | sort -k1n -k2n -k3n -k4n | uniq -c > ${OUT_PREFIX}_gam_dist_h${i}.txt; gzip ${OUT_PREFIX}_gam_dist_h${i}.txt"
done 

# Upload statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/stats/polya_rna/sim_vg_r2/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*txt.gz" --no-progress
