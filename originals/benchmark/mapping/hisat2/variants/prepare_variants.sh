set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"
OUT_PREFIX="${VARIANTS}_${CHR}"

# Download genome
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa.fai . --no-progress

# Download variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${VARIANTS}/${CHR}/${OUT_PREFIX}.vcf.gz variants.vcf.gz --no-progress

# Convert to bi-allelic
/usr/bin/time -v bash -c "bcftools norm -c x -m -any -f ${GENOME_PREFIX}.fa variants.vcf.gz > variants_biallelic.vcf"

if [ "${CHR}" = "X" ] || [ "${CHR}" = "Y" ]; then

	# Convert genotypes to diploid
	/usr/bin/time -v bash -c "head -n 10000 variants_biallelic.vcf | grep '^#' > variants_biallelic_tmp.vcf; grep -v '^#' variants_biallelic.vcf | awk -v OFS='\t' '{for (i = 10; i <= NF; i++) {if (("'$i !~ "\\/"'") && ("'$i !~ "\\|"'")) "'$i'" = "'$i "|" $i'"} {print}}' >> variants_biallelic_tmp.vcf; mv variants_biallelic_tmp.vcf variants_biallelic.vcf"
fi

# Construct variant and haplotype lists
/usr/bin/time -v bash -c "grep -v '^#' variants_biallelic.vcf | wc -l; hisat2_extract_snps_haplotypes_VCF.py --non-rs ${GENOME_PREFIX}.fa variants_biallelic.vcf ${OUT_PREFIX}; wc -l ${OUT_PREFIX}.snp; wc -l ${OUT_PREFIX}.haplotype"

# Upload variants and haplotype lists
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/variants/${VARIANTS}/${CHR} --exclude "*" --include "${OUT_PREFIX}*" --no-progress
