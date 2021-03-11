set -e

# Set genome prefixes
GENOME_SCA_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds"
GENOME_CHR_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"

# Set transcripts prefix
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set variants prefix
VARIANTS_PREFIX="1kg_nonCEU_af001_${CHR}"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_${CHR}"

# Set number of threads
CPU=1

# Chromosome is scaffold
if [ "${CHR}" = "SCA" ]; then

	# Construct variation graph
	/usr/bin/time -v bash -c "vg construct -p -t ${CPU} -r ${GENOME_SCA_PREFIX}.fa > ${CHR}.vg"

	# Convert variation graph to PackedGraph format
	/usr/bin/time -v bash -c "vg convert -p ${CHR}.vg > ${CHR}.pg"

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep '^KI\|^GL' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Construct spliced variation graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -r -n ${CHR}.gtf ${CHR}.pg > ${OUT_PREFIX}.pg"


# Chromosome is mitochondria
elif [ "${CHR}" = "MT" ]; then

	# Construct variation graph
	/usr/bin/time -v bash -c "vg construct -p -t ${CPU} -R ${CHR} -C -r ${GENOME_CHR_PREFIX}.fa > ${CHR}.vg"

	# Convert variation graph to PackedGraph format
	/usr/bin/time -v bash -c "vg convert -p ${CHR}.vg > ${CHR}.pg"

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Construct spliced variation graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -r -n ${CHR}.gtf ${CHR}.pg > ${OUT_PREFIX}.pg"

# Chromosome is canonical
else 

	# Construct variation graph
	/usr/bin/time -v bash -c "vg construct -p -t ${CPU} -R ${CHR} -C -a -v ${VARIANTS_PREFIX}.vcf.gz -r ${GENOME_CHR_PREFIX}.fa > ${CHR}.vg"

	# Convert variation graph to PackedGraph format
	/usr/bin/time -v bash -c "vg convert -p ${CHR}.vg > ${CHR}.pg"

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Construct spliced variation graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -e -n ${CHR}.gtf ${CHR}.pg > ${OUT_PREFIX}.pg"

fi
