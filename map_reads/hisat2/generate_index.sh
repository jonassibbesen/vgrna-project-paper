set -e

# Set genome (FASTA) prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set transcripts (GTF) prefix
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set variants (VCF) prefix
VARIANTS_PREFIX="1kg_nonCEU_af001"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_index"

# Set number of threads
CPU=1

# Construct exon list
/usr/bin/time -v bash -c "hisat2_extract_exons.py ${TRANSCRIPTS_PREFIX}.gtf > ${OUT_PREFIX}_exons.txt; wc -l ${OUT_PREFIX}_exons.txt"

# Construct splice-site list
/usr/bin/time -v bash -c "hisat2_extract_splice_sites.py ${TRANSCRIPTS_PREFIX}.gtf > ${OUT_PREFIX}_splice_sites.txt; wc -l ${OUT_PREFIX}_splice_sites.txt"

if [ "${VARIANTS_PREFIX}" = "" ]; then

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} --exon ${OUT_PREFIX}_exons.txt --ss ${OUT_PREFIX}_splice_sites.txt ${GENOME_PREFIX}.fa ${OUT_PREFIX}"

else

	# Combine variants and haplotype lists
	/usr/bin/time -v bash -c "cat ${VARIANTS_PREFIX}*.snp > ${OUT_PREFIX}_variants.snp; cat ${VARIANTS_PREFIX}*.haplotype > ${OUT_PREFIX}_variants.haplotype; wc -l ${OUT_PREFIX}_variants.snp; wc -l ${OUT_PREFIX}_variants.haplotype"

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} --snp ${OUT_PREFIX}_variants.snp --haplotype ${OUT_PREFIX}_variants.haplotype --exon ${OUT_PREFIX}_exons.txt --ss ${OUT_PREFIX}_splice_sites.txt ${GENOME_PREFIX}.fa ${OUT_PREFIX}"
fi
