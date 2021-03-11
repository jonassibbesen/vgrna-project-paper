set -e

# Set genome prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set output name prefix
OUT_SCA_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds"
OUT_CHR_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"

# Set number of threads
CPU=1

# Select scaffold contigs and index
samtools faidx ${GENOME_PREFIX}.fa $(grep "KI\|GL" ${GENOME_PREFIX}.fa.fai | cut -f1) > ${OUT_SCA_PREFIX}.fa; samtools faidx ${OUT_SCA_PREFIX}.fa; samtools dict ${OUT_SCA_PREFIX}.fa > ${OUT_SCA_PREFIX}.dict

# Select chromosomes + mitochondria and index
samtools faidx ${GENOME_PREFIX}.fa $(grep -v "KI\|GL" ${GENOME_PREFIX}.fa.fai | cut -f1) > ${OUT_CHR_PREFIX}.fa; samtools faidx ${OUT_CHR_PREFIX}.fa; samtools dict ${OUT_CHR_PREFIX}.fa > ${OUT_CHR_PREFIX}.dict
