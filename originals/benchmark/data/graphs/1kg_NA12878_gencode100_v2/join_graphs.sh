set -e

# Set file name prefixes 
GRAPHS_PREFIX="1kg_NA12878_gencode100_v2"

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.pg" --no-progress

# Calculate graph statistics (pre-ids) 
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg stats -z -l -r ${i}/'"${GRAPHS_PREFIX}"'_${i}.pg; done'

# Create non-conflicting id space in graphs
/usr/bin/time -v bash -c 'vg ids -j $(for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}/'"${GRAPHS_PREFIX}"'_${i}.pg; done) -m empty_node_mapping.map'

# Calculate graph statistics (post-ids) 
/usr/bin/time -v bash -c 'for i in $(seq 1 22; echo X; echo Y; echo MT; echo SCA); do echo ${i}; vg stats -z -l -r ${i}/'"${GRAPHS_PREFIX}"'_${i}.pg; done'

# Upload re-id'ed graphs and empty node mapping
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ --exclude "*" --include "*.pg" --no-progress
aws s3 cp empty_node_mapping.map s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ --no-progress
