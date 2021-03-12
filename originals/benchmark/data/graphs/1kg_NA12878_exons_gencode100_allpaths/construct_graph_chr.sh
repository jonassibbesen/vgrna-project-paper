set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"
VARIANTS_PREFIX="1kg_NA12878_exons"
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download genome
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa.fai . --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/${VARIANTS_PREFIX}_${CHR}.vcf.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/${VARIANTS_PREFIX}_${CHR}.vcf.gz.tbi . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Construct variation graph
/usr/bin/time -v bash -c "vg construct -p -t ${CPU} -R ${CHR} -C -a -v ${VARIANTS_PREFIX}_${CHR}.vcf.gz -r ${GENOME_PREFIX}.fa > ${CHR}.vg"

# Convert variation graph to PackedGraph format
/usr/bin/time -v bash -c "vg convert -p ${CHR}.vg > ${CHR}.pg"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Construct spliced variation graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -n ${CHR}.gtf ${CHR}.pg > ${OUT_PREFIX}.pg"

# Upload spliced variation graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
