set -e

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

# Merge gbwts
/usr/bin/time -v bash -c 'vg gbwt -p -m -f -o '"${GRAPH_PREFIX}"'.gbwt $(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/'"${GRAPH_PREFIX}"'_${i}.gbwt; done)'

# Merge seqeunces
/usr/bin/time -v bash -c 'cat $(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/'"${GRAPH_PREFIX}"'_${i}.fa.gz; done) > '"${GRAPH_PREFIX}"'.fa.gz'

# Merge transcript info
/usr/bin/time -v bash -c 'cat $(for i in $(seq 1 22; echo X; echo Y); do echo ${i}/'"${GRAPH_PREFIX}"'_${i}.txt.gz; done) | grep -v ^Name > '"${GRAPH_PREFIX}"'.txt.gz'
