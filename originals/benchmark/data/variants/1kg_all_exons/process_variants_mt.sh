set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"
ALLVARS_PREFIX="ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes"
VARSINFO_PREFIX="integrated_call_samples_v3.20130502.ALL"
EXONS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full_exons"
VARIANTS_PREFIX="1kg_all_exons"
OUT_PREFIX="${VARIANTS_PREFIX}_MT"

# Download genome
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa.fai . --no-progress

# Download variants 
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/variants/1kg_GRCh37_phase3/${ALLVARS_PREFIX}.vcf.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/variants/1kg_GRCh37_phase3/${ALLVARS_PREFIX}.vcf.gz.tbi . --no-progress

# Download info 
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/meta/1000g/${VARSINFO_PREFIX}.panel . --no-progress

# Download exons
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${EXONS_PREFIX}.bed . --no-progress

# Create all samples file
/usr/bin/time -v bash -c "cat ${VARSINFO_PREFIX}.panel | grep -v sample | cut -f1 > samples.txt; wc -l samples.txt"

# Select all samples and exonic variants
/usr/bin/time -v bash -c "bcftools annotate -x ^INFO/AC,^INFO/AF,^INFO/AN ${ALLVARS_PREFIX}.vcf.gz | bcftools norm -m -any -O z -f ${GENOME_PREFIX}.fa > all_variants.vcf.gz; tabix all_variants.vcf.gz; bcftools view -R ${EXONS_PREFIX}.bed -S samples.txt --force-samples -c 1 all_variants.vcf.gz | bcftools norm -m +any -O z -f ${GENOME_PREFIX}.fa > ${OUT_PREFIX}.vcf.gz; tabix ${OUT_PREFIX}.vcf.gz"

# Upload processed variants
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/MT/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
