set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
BASE_PREFIX="1kg_all_af001_gencode100"
GRAPHS_PREFIX="${BASE_PREFIX}_genes"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${CHR}/${BASE_PREFIX}_${CHR}.pg . --no-progress

# Chromosome is scaffold
if [ "${CHR}" = "SCA" ]; then

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep '^KI\|^GL' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${BASE_PREFIX}_${CHR}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -n ${CHR}.gtf ${BASE_PREFIX}_${CHR}.pg > ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Chromosome is mitochondria
elif [ "${CHR}" = "MT" ]; then

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${BASE_PREFIX}_${CHR}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -n ${CHR}.gtf ${BASE_PREFIX}_${CHR}.pg > ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Chromosome is canonical
else 

	HAP_PREFIX="1kg_all_exons"

	# Download haplotypes
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${CHR}/${HAP_PREFIX}_${CHR}.gbwt . --no-progress

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${BASE_PREFIX}_${CHR}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -g -n ${CHR}.gtf -l ${HAP_PREFIX}_${CHR}.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${BASE_PREFIX}_${CHR}.pg > ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Compress haplotype-specific transcripts
	/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

fi

# Upload haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
