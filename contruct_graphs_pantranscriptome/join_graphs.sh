set -e

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

# Calculate graph statistics (pre-ids) 
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg stats -z -l -r ${i}/'"${GRAPH_PREFIX}"'_${i}.pg; done'

# Create non-conflicting id space in graphs
/usr/bin/time -v bash -c 'vg ids -j $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}/'"${GRAPH_PREFIX}"'_${i}.pg; done) -m empty_node_mapping.map'

# Calculate graph statistics (post-ids) 
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg stats -z -l -r ${i}/'"${GRAPH_PREFIX}"'_${i}.pg; done'
