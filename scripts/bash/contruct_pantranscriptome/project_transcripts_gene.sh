set -e

# Set transcripts (GTF) prefix
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"

# Set variants (VCF) prefix
VARIANTS_PREFIX="1kg_nonCEU_exons_${CHR}"

# Set graph (PG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100_${CHR}"

# Set output name prefix
OUT_PREFIX="1kg_nonCEU_af001_gencode100_genes_${CHR}"

# Set number of threads
CPU=1

# Chromosome is scaffold
if [ "${CHR}" = "SCA" ]; then

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep '^KI\|^GL' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${GRAPH_PREFIX}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -n ${CHR}.gtf ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Chromosome is mitochondria
elif [ "${CHR}" = "MT" ]; then

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${GRAPH_PREFIX}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -n ${CHR}.gtf ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}.pg"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Chromosome is canonical
else 

	# Create gbwt index of all haplotypes
	/usr/bin/time -v bash -c "vg index -p -t ${CPU} -G haplotypes.gbwt -v ${VARIANTS_PREFIX}.vcf.gz ${GRAPH_PREFIX}.pg"

	# Find contig transcripts
	/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

	# Calculate graph statistics (pre-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${GRAPH_PREFIX}.pg"

	# Create haplotype-specific transcripts and update graph
	/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -d -o -r -g -n ${CHR}.gtf -l haplotypes.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${GRAPH_PREFIX}.pg > ${OUT_PREFIX}.pg; rm haplotypes.gbwt"

	# Calculate graph statistics (post-rna) 
	/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

	# Compress haplotype-specific transcripts
	/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

fi

# Upload haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
