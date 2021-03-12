set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
BASE_PREFIX="1kg_EURnonCEU_af002_gencode100"
HAP_PREFIX="1kg_EURnonCEU_exons"
GRAPHS_PREFIX="${BASE_PREFIX}_unidi"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${CHR}/${BASE_PREFIX}_${CHR}.pg . --no-progress

# Download haplotypes
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${CHR}/${HAP_PREFIX}_${CHR}.gbwt . --no-progress

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${BASE_PREFIX}_${CHR}.pg"

# Create haplotype-specific transcripts and update graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -n ${CHR}.gtf -l ${HAP_PREFIX}_${CHR}.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${BASE_PREFIX}_${CHR}.pg > ${BASE_PREFIX}_${CHR}_2.pg"

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${BASE_PREFIX}_${CHR}_2.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Upload haplotype-specific transcripts
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
