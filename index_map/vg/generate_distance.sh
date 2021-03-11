set -e

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_index"

# Set number of threads
CPU=1

# Generate trivial snarls
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg snarls -t '"${CPU}"' --algorithm integrated -T ${i}/'"${GRAPH_PREFIX}"'_${i}.pg > ${i}_trivial.snarls; done'

# Combine trivial snarls
/usr/bin/time -v bash -c 'cat $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}_trivial.snarls; done) > trivial.snarls; rm *_trivial.snarls'

# Generate distance index
/usr/bin/time -v bash -c "vg index -p -t ${CPU} -x ${GRAPH_PREFIX}.xg -s trivial.snarls -j ${OUT_PREFIX}.dist"
