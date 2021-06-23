set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"
VARIANTS_PREFIX="1kg_all_exons"
OUT_PREFIX="${VARIANTS_PREFIX}_${CHR}_vep"

# Download VEP cache
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/databases/vep/homo_sapiens_vep_103_GRCh38.tar.gz . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download genome
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa.fai . --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/${VARIANTS_PREFIX}_${CHR}.vcf.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/${VARIANTS_PREFIX}_${CHR}.vcf.gz.tbi . --no-progress

# Unpack VEP cache 
/usr/bin/time -v bash -c "tar -xzf homo_sapiens_vep_103_GRCh38.tar.gz"

# Find contig transcripts, sort and index
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ${CHR}.gtf.gz; tabix -p gff ${CHR}.gtf.gz"

# Find chromosome seqeunce and index
/usr/bin/time -v bash -c "samtools faidx ${GENOME_PREFIX}.fa ${CHR} > ${CHR}.fa; samtools faidx ${CHR}.fa"

# Predict variant consequences
/usr/bin/time -v bash -c "/opt/vep/src/ensembl-vep/vep --verbose --species homo_sapiens --everything --check_existing --offline --tab --minimal --allele_number --show_ref_allele --cache --dir_cache . -i ${VARIANTS_PREFIX}_${CHR}.vcf.gz --gtf ${CHR}.gtf.gz --fasta ${CHR}.fa -o ${OUT_PREFIX}.txt"

# Compress predictions
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt"

# Upload processed variants
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
