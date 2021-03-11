set -e

# Set genome prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set variants prefix
VARIANTS_PREFIX="1kg_nonCEU_af001_${CHR}"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_${CHR}"

# Set number of threads
CPU=1

# Convert to bi-allelic
/usr/bin/time -v bash -c "bcftools norm -c x -m -any -f ${GENOME_PREFIX}.fa ${VARIANTS_PREFIX}.vcf.gz > variants_biallelic.vcf"

if [ "${CHR}" = "X" ] || [ "${CHR}" = "Y" ]; then

	# Convert genotypes to diploid
	/usr/bin/time -v bash -c "head -n 10000 variants_biallelic.vcf | grep '^#' > variants_biallelic_tmp.vcf; grep -v '^#' variants_biallelic.vcf | awk -v OFS='\t' '{for (i = 10; i <= NF; i++) {if (("'$i !~ "\\/"'") && ("'$i !~ "\\|"'")) "'$i'" = "'$i "|" $i'"} {print}}' >> variants_biallelic_tmp.vcf; mv variants_biallelic_tmp.vcf variants_biallelic.vcf"
fi

# Construct variant and haplotype lists
/usr/bin/time -v bash -c "grep -v '^#' variants_biallelic.vcf | wc -l; hisat2_extract_snps_haplotypes_VCF.py --non-rs ${GENOME_PREFIX}.fa variants_biallelic.vcf ${OUT_PREFIX}; wc -l ${OUT_PREFIX}.snp; wc -l ${OUT_PREFIX}.haplotype; rm variants_biallelic.vcf"
