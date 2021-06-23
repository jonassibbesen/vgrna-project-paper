set -e

# Set file name prefixes
ALIGN_PREFIX="mpmap_${GRAPH}_real_${NAME}"
OUT_PREFIX="rpvg_mpmap_strand_${NAME}_${TRANSCRIPTS}"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/alignments/polya_rna/real/${NAME}/mpmap/${GRAPH}/${ALIGN_PREFIX}.gamp . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt . --no-progress

# Download transcript info
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Concatenate transcript info
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep -v ^Name > transcript_info.txt; wc -l transcript_info.txt"

# Infer expression
/usr/bin/time -v bash -c "/rpvg/bin/rpvg -t ${CPU} -r ${SEED} -e rf -n 1 -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a ${ALIGN_PREFIX}.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt; gzip ${OUT_PREFIX}_haps.txt"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/analyses/imprinting/ENCSR000AED/inference/${NAME}/${TRANSCRIPTS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
