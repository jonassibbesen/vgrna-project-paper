set -e

MAF_SUFFIX=$(echo ${MAF} | cut -d '.' -f2)

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"
ALLVARS_PREFIX="ALL.chr${CHR}_GRCh38.genotypes.20170504"
VARSINFO_PREFIX="integrated_call_samples_v3.20130502.ALL"
EXONVARS_PREFIX="1kg_EURnonCEU_exons"
VARIANTS_PREFIX="1kg_EURnonCEU_af${MAF_SUFFIX}"
OUT_PREFIX="${VARIANTS_PREFIX}_${CHR}"

# Download genome
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa.fai . --no-progress

# Download variants 
aws s3 cp s3://vg-data/1kg_GRCh38/variants/${ALLVARS_PREFIX}.vcf.gz . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/variants/${ALLVARS_PREFIX}.vcf.gz.tbi . --no-progress

# Download info 
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/meta/1000g/${VARSINFO_PREFIX}.panel . --no-progress

# Download exon variants 
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${EXONVARS_PREFIX}/${CHR}/${EXONVARS_PREFIX}_${CHR}.vcf.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${EXONVARS_PREFIX}/${CHR}/${EXONVARS_PREFIX}_${CHR}.vcf.gz.tbi . --no-progress

# Create EUR-non-CEU samples file
/usr/bin/time -v bash -c "grep EUR ${VARSINFO_PREFIX}.panel | grep -v CEU | grep -v sample | cut -f1 > samples.txt; wc -l samples.txt"

# Add female samples to chromosome Y 
if [ "${CHR}" = "Y" ]; then

	aws s3 cp s3://vg-data/1kg_GRCh38/variants/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz .

	/usr/bin/time -v bash -c "zcat ALL.chr22_GRCh38.genotypes.20170504.vcf.gz | head -n 1000 | grep ^# > header.vcf; bgzip header.vcf; tabix header.vcf.gz; bcftools merge --force-samples -O z ${ALLVARS_PREFIX}.vcf.gz header.vcf.gz > ${ALLVARS_PREFIX}_tmp.vcf.gz; mv ${ALLVARS_PREFIX}_tmp.vcf.gz ${ALLVARS_PREFIX}.vcf.gz; tabix ${ALLVARS_PREFIX}.vcf.gz"
fi

# Select EUR-non-CEU samples
/usr/bin/time -v bash -c "bcftools annotate -x ^INFO/AC,^INFO/AF,^INFO/AN ${ALLVARS_PREFIX}.vcf.gz | bcftools norm -m -any -f ${GENOME_PREFIX}.fa | bcftools view -S samples.txt -q ${MAF} | bcftools norm -m +any -O z -f ${GENOME_PREFIX}.fa > all_variants.vcf.gz; tabix all_variants.vcf.gz"

# Concat variant sets
/usr/bin/time -v bash -c "bcftools concat -a -d all -O z ${EXONVARS_PREFIX}_${CHR}.vcf.gz all_variants.vcf.gz > ${OUT_PREFIX}.vcf.gz; tabix ${OUT_PREFIX}.vcf.gz"

# Upload processed variants
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
