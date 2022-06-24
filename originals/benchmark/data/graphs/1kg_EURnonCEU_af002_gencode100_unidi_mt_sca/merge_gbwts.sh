set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_EURnonCEU_af002_gencode100_unidi_mt_sca"
OUT_PREFIX="${GRAPHS_PREFIX}"

# Download gbwts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.gbwt" --no-progress

# Merge gbwts
/usr/bin/time -v bash -c 'vg gbwt -p -m -f -o '"${GRAPHS_PREFIX}"'.gbwt $(for i in $(echo MT; echo SCA); do echo ${i}/'"${GRAPHS_PREFIX}"'_${i}.gbwt; done)'

# Upload merged gbwt
aws s3 cp ${OUT_PREFIX}.gbwt s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ --no-progress
