set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_nonCEU_af001_imgt_hla_p10k_noB258_noN_main_gencode100"
OUT_PREFIX="${GRAPHS_PREFIX}"

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.pg" --no-progress

# Generate xg graph
/usr/bin/time -v bash -c 'vg index -p -t '"${CPU}"' -x '"${GRAPHS_PREFIX}"'.xg $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}/'"${GRAPHS_PREFIX}"'_${i}.pg; done); vg paths -L -x '"${OUT_PREFIX}"'.xg | wc -l'

# Upload xg graph
aws s3 cp ${OUT_PREFIX}.xg s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ --no-progress
