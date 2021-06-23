set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
GRAPHS_PREFIX="1kg_nonCEU_af001_imgt_hla_p10k_noB258_noN_main_gencode100"
OUT_PREFIX="${GRAPHS_PREFIX}_6"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf transcripts.gtf --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/${GRAPHS_PREFIX}_6.pg 6.pg --no-progress

# Extract haplotype path
/usr/bin/time -v bash -c "vg paths -Q hla -v 6.pg -A > haps.gaf; wc -l haps.gaf"

# Remove haplotype path
/usr/bin/time -v bash -c "vg paths -L -v 6.pg | wc -l; vg paths -Q hla -v 6.pg -d | vg convert -t ${CPU} -p - > 6_nohaps.pg; vg paths -L -v 6_nohaps.pg | wc -l"

# Convert haplotype paths to GBWT
/usr/bin/time -v bash -c "vg gbwt -p --num-threads ${CPU} -A haps.gaf -x 6_nohaps.pg -o ${OUT_PREFIX}_haps.gbwt"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^6\t' transcripts.gtf > 6.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r 6_nohaps.pg"

# Create haplotype-specific transcripts and update graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -g -n 6.gtf -l ${OUT_PREFIX}_haps.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt 6_nohaps.pg > ${OUT_PREFIX}.pg"

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Upload haplotypes, haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
