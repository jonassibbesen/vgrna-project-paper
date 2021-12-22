set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="1kg_NA12878_gencode100"
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
READS_PREFIX="SRR1153470"
EXPRESSION_PREFIX="${READS_PREFIX}_uni"
PROFILE_PREFIX="${TRANSCRIPTS_PREFIX}_${EXPRESSION_PREFIX}"
MODEL_PREFIX="${TRANSCRIPTS_PREFIX}_${READS_PREFIX}_rsem"
OUT_PREFIX="sim_${PROFILE_PREFIX}_rsem_r1_h1"

# Download haplotype-specific transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.fa.gz" --include "*.txt.gz" --no-progress

# Download expression profile
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/${TRANSCRIPTS_PREFIX}/${EXPRESSION_PREFIX}/${PROFILE_PREFIX}.isoforms.results . --no-progress

# Download simulation model
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/${TRANSCRIPTS_PREFIX}/${READS_PREFIX}/${MODEL_PREFIX}.stat/${MODEL_PREFIX}.model . --no-progress

# Extract haplotype-specific transcripts
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep NA12878 | grep _0_ | cut -f1 | grep -v Name > haplotype_transcripts.txt; seqtk subseq <(zcat */*fa.gz) haplotype_transcripts.txt > ${TRANSCRIPTS_PREFIX}.fa; grep '"'>'"' ${TRANSCRIPTS_PREFIX}.fa | wc -l"

# Generate transcriptome reference
/usr/bin/time -v bash -c "rsem-prepare-reference -p ${CPU} ${TRANSCRIPTS_PREFIX}.fa ${TRANSCRIPTS_PREFIX}"

# Extract haplotype-specific expression profile
/usr/bin/time -v bash -c "echo transcript_id >> haplotype_transcripts.txt; grep -F -f haplotype_transcripts.txt ${PROFILE_PREFIX}.isoforms.results > ${PROFILE_PREFIX}_haplotype.isoforms.results; wc -l ${PROFILE_PREFIX}_haplotype.isoforms.results"

# Simulate reads
/usr/bin/time -v bash -c "rsem-simulate-reads ${TRANSCRIPTS_PREFIX} ${MODEL_PREFIX}.model ${PROFILE_PREFIX}_haplotype.isoforms.results 0 ${NREADS} ${OUT_PREFIX} --seed ${SEED}; wc -l ${OUT_PREFIX}_1.fq; wc -l ${OUT_PREFIX}_2.fq"

# Compress simulated reads
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}_1.fq; gzip ${OUT_PREFIX}_2.fq"

# Upload simulated reads
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/${TRANSCRIPTS_PREFIX}/${EXPRESSION_PREFIX}/rsem_r1/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
