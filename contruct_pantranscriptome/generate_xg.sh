set -e

# Set graph (PG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

# Generate xg graph
/usr/bin/time -v bash -c 'vg index -p -t '"${CPU}"' -x '"${GRAPH_PREFIX}"'.xg $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo '"${GRAPH_PREFIX}"'_${i}.pg; done); vg paths -L -x '"${OUT_PREFIX}"'.xg | wc -l'
