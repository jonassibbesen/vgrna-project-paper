set -e

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_index"

# Set number of threads
CPU=1

# Prune graphs
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg prune -p -t '"${CPU}"' -r ${i}/'"${GRAPH_PREFIX}"'_${i}.pg > '"${GRAPH_PREFIX}"'_${i}_pruned.vg; done'

# Generate gcsa index
/usr/bin/time -v bash -c 'vg index -p -t '"${CPU}"' -g '"${OUT_PREFIX}.gcsa"' $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo '"${GRAPH_PREFIX}"'_${i}_pruned.vg; done); rm *_pruned.vg'
