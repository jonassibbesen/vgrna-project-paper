set -e

# Set file name prefixes
READS_PREFIX="sim_1kg_NA12878_gencode100_${SIM}_rsem_r1"
OUT_PREFIX="${MAPPER}_${REF}_sim_rsem_r1_${SIM}"

# Download reads
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/${SIM}/rsem_r1/ . --recursive --exclude "*" --include "*.fq.gz" --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${REF}/ . --recursive --exclude "*" --include "*.xg" --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${REF}/ . --recursive --no-progress

for i in $(seq 1 2); do

	# Use default vg map
	if [ "${MAPPER}" = "map" ]; then

		# Map reads
		/usr/bin/time -v bash -c "vg map -t ${CPU} -x ${REF}.xg -g ${REF}_index.gcsa -f ${READS_PREFIX}_h${i}_1.fq.gz -f ${READS_PREFIX}_h${i}_2.fq.gz > ${OUT_PREFIX}_h${i}.gam"

	# Use faster vg map
	elif [ "${MAPPER}" = "map_fast" ]; then

		# Map reads
		/usr/bin/time -v bash -c "vg map -t ${CPU} --try-up-to 16 --mate-rescues 32 -x ${REF}.xg -g ${REF}_index.gcsa -f ${READS_PREFIX}_h${i}_1.fq.gz -f ${READS_PREFIX}_h${i}_2.fq.gz > ${OUT_PREFIX}_h${i}.gam"

	# Use default vg mpmap
	elif [ "${MAPPER}" = "mpmap" ]; then

		# Map reads
		/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna --report-allelic-mapq -x ${REF}.xg -g ${REF}_index.gcsa -d ${REF}_index.dist -f ${READS_PREFIX}_h${i}_1.fq.gz -f ${READS_PREFIX}_h${i}_2.fq.gz > ${OUT_PREFIX}_h${i}.gamp"

	# Use vg mpmap in non-splicing mode
	elif [ "${MAPPER}" = "mpmap_nosplice" ]; then

		# Map reads
		/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna --report-allelic-mapq --not-spliced -x ${REF}.xg -g ${REF}_index.gcsa -d ${REF}_index.dist -f ${READS_PREFIX}_h${i}_1.fq.gz -f ${READS_PREFIX}_h${i}_2.fq.gz > ${OUT_PREFIX}_h${i}.gamp"
	fi
done

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/sim_rsem_r1/${SIM}/${MAPPER}/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
