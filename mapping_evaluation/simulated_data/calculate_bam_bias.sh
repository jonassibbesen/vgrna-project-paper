set -e

# Set alignment prefixes
ALIGN1_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"
ALIGN2_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h2"

# Set variants prefix
VARIANTS_PREFIX="1kg_NA12878_exons_${CHR}"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_${CHR}"

# Set number of threads
CPU=1

# Calculate allele coverage
/usr/bin/time -v bash -c "calc_allele_read_coverage ${ALIGN1_PREFIX}.bam ${ALIGN2_PREFIX}.bam ${VARIANTS_PREFIX}.vcf > ${OUT_PREFIX}_allele_cov.txt; gzip ${OUT_PREFIX}_allele_cov.txt"	 
