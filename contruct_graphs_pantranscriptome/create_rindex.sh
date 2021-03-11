set -e

# Set input name prefixes
GRAPH_PREFIX="1kg_all_af001_gencode100"

# Set output name prefixes
OUT_PREFIX="1kg_all_af001_gencode100"

# Set number of threads
CPU=1

# Create r-index
/usr/bin/time -v bash -c "vg gbwt --num-threads ${CPU} -r ${OUT_PREFIX}.gbwt.ri ${GRAPH_PREFIX}.gbwt"
