set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
GRAPHS_PREFIX="gencode100"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/${OUT_PREFIX}.pg . --no-progress

# Chromosome is scaffold
if [ "${CHR}" = "SCA" ]; then

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep '^KI\|^GL' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -o -r -u -g -n ${CHR}.gtf -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${OUT_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Compress haplotype-specific transcripts
	/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Chromosome is mitochondria or canonical
else 

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -o -r -u -g -n ${CHR}.gtf -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${OUT_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Compress haplotype-specific transcripts
	/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

fi

# Upload haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
