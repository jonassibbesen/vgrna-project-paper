set -e

# Set file name prefixes 
OUT_PREFIX="${GRAPHS}_index"

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/ . --recursive --exclude "*" --include "*.pg" --no-progress

# Prune graphs
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg prune -p -t '"${CPU}"' -r ${i}/'"${GRAPHS}"'_${i}.pg > '"${GRAPHS}"'_${i}_pruned.vg; done'

# Generate gcsa index
/usr/bin/time -v bash -c 'vg index -p -t '"${CPU}"' -g '"${OUT_PREFIX}.gcsa"' $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo '"${GRAPHS}"'_${i}_pruned.vg; done)'

# Upload gcsa index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${GRAPHS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
