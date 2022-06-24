set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
BASE_PREFIX="1kg_all_af001_gencode100"
GRAPHS_PREFIX="${BASE_PREFIX}_unidi_mt_sca"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${CHR}/${BASE_PREFIX}_${CHR}.pg graph.pg --no-progress

# Chromosome is scaffold
if [ "${CHR}" = "SCA" ]; then

    # Find contig transcripts
    /usr/bin/time -v bash -c "grep '^KI\|^GL' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Chromosome is mitochondria
else 

	# Find contig transcripts
    /usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"
fi

# Remove transcript paths
/usr/bin/time -v bash -c "vg paths -d -Q ENST -v graph.pg > graph_nopaths.pg"

# Create reference transcripts
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -o -u -n ${CHR}.gtf -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt graph_nopaths.pg > graph_nopaths2.pg"

# Compress reference transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Upload reference transcripts
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
