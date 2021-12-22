set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_all_af001_mt_gencode100_unidi"
OUT_PREFIX="${GRAPHS_PREFIX}_main"

# Download gbwts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${GRAPHS_PREFIX}.gbwt . --no-progress

# Create r-index
/usr/bin/time -v bash -c "vg gbwt --num-threads ${CPU} -r ${OUT_PREFIX}.gbwt.ri ${GRAPHS_PREFIX}.gbwt"

# Upload r-index
aws s3 cp ${OUT_PREFIX}.gbwt.ri s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${OUT_PREFIX}/ --no-progress
