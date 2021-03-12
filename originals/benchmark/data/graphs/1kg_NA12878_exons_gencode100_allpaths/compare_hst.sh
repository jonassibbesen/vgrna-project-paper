set -e

# Set file name prefixes 
HST_1="1kg_NA12878_exons_gencode100_allpaths"
OUT_PREFIX="${HST_2}_hst_overlap"

# Download hst sequnces 1
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${HST_1}/ . --recursive --exclude "*" --include "*.fa.gz" --no-progress

# Download hst sequnces 2
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${HST_2}/ . --recursive --exclude "*" --include "*.fa.gz" --no-progress

# Concatenate hst sequences 1
/usr/bin/time -v bash -c "zcat */${HST_1}*fa.gz > ${HST_1}.fa; cat ${HST_1}.fa | grep '"'>'"' | wc -l"

# Concatenate hst sequences 2
/usr/bin/time -v bash -c "zcat */${HST_2}*fa.gz > ${HST_2}.fa; cat ${HST_2}.fa | grep '"'>'"' | wc -l"

# Compare hst sequences
/usr/bin/time -v bash -c "python3 /scripts/vgrna/compare_hst_sequences.py ${HST_1}.fa ${HST_2}.fa ${OUT_PREFIX}.txt; wc -l ${OUT_PREFIX}.txt"

# Upload hst overlap
aws s3 cp ${OUT_PREFIX}.txt s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${HST_1}/ --no-progress
