set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_all_af001_imgt_hla_p10k_noB258_noN_a100_gencode100"
OUT_PREFIX="${GRAPHS_PREFIX}_6"

# Download gbwts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/${GRAPHS_PREFIX}_6.gbwt . --no-progress

# Create r-index
/usr/bin/time -v bash -c "vg gbwt --num-threads ${CPU} -r ${OUT_PREFIX}.gbwt.ri ${GRAPHS_PREFIX}_6.gbwt"

# Upload r-index
aws s3 cp ${OUT_PREFIX}.gbwt.ri s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/ --no-progress
