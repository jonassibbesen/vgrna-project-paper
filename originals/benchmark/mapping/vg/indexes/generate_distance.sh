set -e

# Set file name prefixes 
OUT_PREFIX="${GRAPHS}_index"

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/ . --recursive --exclude "*" --include "*.pg" --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS}/${GRAPHS}.xg . --no-progress

# Generate trivial snarls
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg snarls -t '"${CPU}"' --algorithm integrated -T ${i}/'"${GRAPHS}"'_${i}.pg > ${i}_trivial.snarls; done'

# Combine trivial snarls
/usr/bin/time -v bash -c 'cat $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}_trivial.snarls; done) > trivial.snarls'

# Generate distance index
/usr/bin/time -v bash -c "vg index -p -t ${CPU} -x ${GRAPHS}.xg -s trivial.snarls -j ${OUT_PREFIX}.dist"

# Upload distance index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${GRAPHS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
