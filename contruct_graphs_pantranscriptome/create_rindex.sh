set -e

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

# Create r-index
/usr/bin/time -v bash -c "vg gbwt --num-threads ${CPU} -r ${OUT_PREFIX}.gbwt.ri ${GRAPH_PREFIX}.gbwt"
