set -e

# Set file name prefixes 
EXONVARS_PREFIX="1kg_NA12878_exons"
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
GRAPHS_PREFIX="1kg_NA12878_gencode100_v2"
OUT_PREFIX="${GRAPHS_PREFIX}_${CHR}"

# Download exon variants
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${EXONVARS_PREFIX}/${CHR}/${EXONVARS_PREFIX}_${CHR}.vcf.gz . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/${EXONVARS_PREFIX}/${CHR}/${EXONVARS_PREFIX}_${CHR}.vcf.gz.tbi . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf . --no-progress

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/${OUT_PREFIX}.pg . --no-progress

# Create gbwt index of NA12878 haplotypes
/usr/bin/time -v bash -c "vg gbwt -p --num-threads ${CPU} --discard-overlaps -v ${EXONVARS_PREFIX}_${CHR}.vcf.gz -x ${OUT_PREFIX}.pg -o ${EXONVARS_PREFIX}_${CHR}.gbwt"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^${CHR}\t' ${TRANSCRIPTS_PREFIX}.gtf > ${CHR}.gtf"

# Calculate graph statistics (pre-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Create haplotype-specific transcripts and update graph
/usr/bin/time -v bash -c "vg rna -p -t ${CPU} -o -r -g -n ${CHR}.gtf -l ${EXONVARS_PREFIX}_${CHR}.gbwt -b ${OUT_PREFIX}.gbwt -f ${OUT_PREFIX}.fa -i ${OUT_PREFIX}.txt ${OUT_PREFIX}.pg > ${OUT_PREFIX}_tmp.pg; mv ${OUT_PREFIX}_tmp.pg ${OUT_PREFIX}.pg"

# Calculate graph statistics (post-rna) 
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Compress haplotype-specific transcripts
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.fa; gzip ${OUT_PREFIX}.txt"

# Upload haplotypes, haplotype-specific transcripts and updated graph
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${CHR}/ --exclude "*" --include "${EXONVARS_PREFIX}_${CHR}.gbwt" --include "${OUT_PREFIX}*" --no-progress
