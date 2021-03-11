set -e

# Set graph prefix
GRAPH_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set output name prefix
OUT_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set number of threads
CPU=1

# Extract transcripts as alignments
/usr/bin/time -v bash -c "vg paths -L -x ${GRAPH_PREFIX}.xg | wc -l; vg paths -X -Q ENST -x ${GRAPH_PREFIX}.xg > ${GRAPH_PREFIX}.gam; vg view -a ${GRAPH_PREFIX}.gam | wc -l"

# Extract reference path list
/usr/bin/time -v bash -c "vg paths -L -x ${GRAPH_PREFIX}.xg | grep -v ENST > reference_paths.txt; wc -l reference_paths.txt"

# Surject transcript alignments to reference
/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -b -F reference_paths.txt -x ${GRAPH_PREFIX}.xg ${GRAPH_PREFIX}.gam | samtools sort --threads ${CPU} - > ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam; samtools flagstat ${OUT_PREFIX}.bam; rm reference_paths.txt"
