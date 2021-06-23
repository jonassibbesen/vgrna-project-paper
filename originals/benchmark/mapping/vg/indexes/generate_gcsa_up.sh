set -e

# Set file name prefixes 
OUT_PREFIX="${GRAPHS}_index_up"

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/ . --recursive --exclude "*" --include "*.pg" --no-progress

# Download pantranscriptomes
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/ . --recursive --exclude "*" --include "*.gbwt" --no-progress

# Download empty node-mapping
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/empty_node_mapping.map empty_node_mapping.map --no-progress

# Prune graphs
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y); do echo ${i}; vg prune -p -t '"${CPU}"' -u -a -m empty_node_mapping.map -g ${i}/'"${GRAPHS}"'_${i}.gbwt ${i}/'"${GRAPHS}"'_${i}.pg > '"${GRAPHS}"'_${i}_pruned.vg; done'

# Prune graphs
/usr/bin/time -v bash -c 'for i in $(echo MT; echo SCA); do echo ${i}; vg prune -p -t '"${CPU}"' -r ${i}/'"${GRAPHS}"'_${i}.pg > '"${GRAPHS}"'_${i}_pruned.vg; done'

# Generate gcsa index
/usr/bin/time -v bash -c 'vg index -p -t '"${CPU}"' -f empty_node_mapping.map -g '"${OUT_PREFIX}.gcsa"' $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo '"${GRAPHS}"'_${i}_pruned.vg; done)'

# Upload gcsa index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${GRAPHS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
