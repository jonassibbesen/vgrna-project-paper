set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
BASE_PREFIX="1kg_all_af001_imgt_hla_p10k_noB258_noN_a100_gencode100"
GRAPHS_PREFIX="${BASE_PREFIX}_unidi"
OUT_PREFIX="${GRAPHS_PREFIX}_6"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf transcripts.gtf --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/6/${BASE_PREFIX}_6.pg 6.pg --no-progress

# Download haplotypes
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/6/${BASE_PREFIX}_6_haps.gbwt haps.gbwt --no-progress

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^6\t' transcripts.gtf > 6.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r 6.pg"

# Create haplotype-specific transcripts and update graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -n 6.gtf -l haps.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt 6.pg > 6_2.pg"

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r 6_2.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Upload haplotypes, haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
