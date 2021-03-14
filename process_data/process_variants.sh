set -e

# Set minimum allele frequency
MAF="0.001"
MAF_SUFFIX=$(echo ${MAF} | cut -d '.' -f2)

# Set samples file
SAMPLES="samples.txt"

# Set genome (FASTA) prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"

# Set exons (BED) prefix
EXONS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full_exons"

# Set variants (VCF) prefix
VARIANTS_22_PREFIX="ALL.chr22_GRCh38.genotypes.20170504"
VARIANTS_CHR_PREFIX="ALL.chr${CHR}_GRCh38.genotypes.20170504"

# Set output name prefixes
OUT_EXON_PREFIX="1kg_all_exons_${CHR}"
OUT_ALL_PREFIX="1kg_all_af{MAF_SUFFIX}_${CHR}"

# Set number of threads
CPU=1

# Add female samples to chromosome Y 
if [ "${CHR}" = "Y" ]; then

	VARIANTS_CHR_PREFIX="${VARIANTS_CHR_PREFIX}_female"

	zcat ${VARIANTS_22_PREFIX}.vcf.gz | head -n 1000 | grep ^# > header.vcf; bgzip header.vcf; tabix header.vcf.gz; bcftools merge --force-samples -O z ${VARIANTS_CHR_PREFIX}.vcf.gz header.vcf.gz > ${VARIANTS_CHR_PREFIX}_tmp.vcf.gz; mv ${VARIANTS_CHR_PREFIX}_tmp.vcf.gz ${VARIANTS_CHR_PREFIX}.vcf.gz; tabix ${VARIANTS_CHR_PREFIX}.vcf.gz
fi

# Select samples
bcftools annotate -x ^INFO/AC,^INFO/AF,^INFO/AN ${VARIANTS_CHR_PREFIX}.vcf.gz | bcftools norm -m -any -f ${GENOME_PREFIX}.fa | bcftools view -S ${SAMPLES} -q ${MAF} | bcftools norm -m +any -O z -f ${GENOME_PREFIX}.fa > all_variants.vcf.gz; tabix all_variants.vcf.gz

# Select samples and exonic variants
bcftools annotate -x ^INFO/AC,^INFO/AF,^INFO/AN ${VARIANTS_CHR_PREFIX}.vcf.gz | bcftools norm -m -any -O z -f ${GENOME_PREFIX}.fa > all_variants.vcf.gz; tabix all_variants.vcf.gz; bcftools view -R ${EXONS_PREFIX}.bed -S ${SAMPLES} -c 1 all_variants.vcf.gz | bcftools norm -m +any -O z -f ${GENOME_PREFIX}.fa > ${OUT_EXON_PREFIX}.vcf.gz; tabix ${OUT_EXON_PREFIX}.vcf.gz

# Concat variant sets
bcftools concat -a -d all -O z ${OUT_EXON_PREFIX}.vcf.gz all_variants.vcf.gz > ${OUT_ALL_PREFIX}.vcf.gz; tabix ${OUT_ALL_PREFIX}.vcf.gz; rm all_variants.vcf.gz*
