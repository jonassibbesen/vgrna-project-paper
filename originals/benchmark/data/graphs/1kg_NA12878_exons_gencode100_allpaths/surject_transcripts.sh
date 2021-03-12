set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
OUT_PREFIX="${GRAPHS_PREFIX}"

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${GRAPHS_PREFIX}.xg . --no-progress

# Extract transcripts as alignments
/usr/bin/time -v bash -c "vg paths -L -x ${GRAPHS_PREFIX}.xg | wc -l; vg paths -X -Q ENST -x ${GRAPHS_PREFIX}.xg > ${GRAPHS_PREFIX}.gam; vg view -a ${GRAPHS_PREFIX}.gam | wc -l"

# Extract reference path list
/usr/bin/time -v bash -c "vg paths -L -x ${GRAPHS_PREFIX}.xg | grep -v ENST > reference_paths.txt; wc -l reference_paths.txt"

# Surject transcript alignments to reference
/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -b -F reference_paths.txt -x ${GRAPHS_PREFIX}.xg ${GRAPHS_PREFIX}.gam | samtools sort --threads ${CPU} - > ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam; samtools flagstat ${OUT_PREFIX}.bam"

# Upload surjected transcripts
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ --exclude "*" --include "${OUT_PREFIX}.bam*" --no-progress
