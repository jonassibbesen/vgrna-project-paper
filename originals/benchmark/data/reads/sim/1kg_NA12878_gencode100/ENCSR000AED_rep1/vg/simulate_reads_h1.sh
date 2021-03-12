set -e

# Set file name prefixes 
TRANSCRIPTS_PREFIX="1kg_NA12878_gencode100"
GRAPHS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"
READS_PREFIX="ENCSR000AED_rep1"
EXPRESSION_PREFIX="${READS_PREFIX}"
PROFILE_PREFIX="${TRANSCRIPTS_PREFIX}_${EXPRESSION_PREFIX}_rsem"
OUT_PREFIX="sim_${TRANSCRIPTS_PREFIX}_${EXPRESSION_PREFIX}_vg_h1"

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/${GRAPHS_PREFIX}.xg . --no-progress

# Download haplotype-specific transcript summary
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Download reads
aws s3 cp s3://encode-public/2013/06/13/c653a32e-e618-42b1-b8b8-b3b838847b97/ENCFF001REK.fastq.gz . --no-sign-request --no-progress
aws s3 cp s3://encode-public/2013/06/13/efa1a02d-6b43-4635-9ef8-d2d78c527839/ENCFF001REJ.fastq.gz . --no-sign-request --no-progress

# Download expression profile
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/${TRANSCRIPTS_PREFIX}/${READS_PREFIX}/${PROFILE_PREFIX}.isoforms.results . --no-progress

# Extract haplotype-specific expression
/usr/bin/time -v bash -c "zcat */*.txt.gz | grep NA12878 | grep _0_ | cut -f1 | grep -v Name > haplotype_transcripts.txt; echo transcript_id >> haplotype_transcripts.txt; wc -l haplotype_transcripts.txt; grep -F -f haplotype_transcripts.txt ${PROFILE_PREFIX}.isoforms.results > ${PROFILE_PREFIX}_haplotype.isoforms.results; wc -l ${PROFILE_PREFIX}_haplotype.isoforms.results"

# Subset read pairs
/usr/bin/time -v bash -c "seqtk sample -s ${SEED} ENCFF001REK.fastq.gz 10000000 > ${READS_PREFIX}_subset_1.fq; seqtk sample -s ${SEED} ENCFF001REJ.fastq.gz 10000000 > ${READS_PREFIX}_subset_2.fq; wc -l ${READS_PREFIX}_subset_1.fq; wc -l ${READS_PREFIX}_subset_2.fq"

# Simulate reads
/usr/bin/time -v bash -c "vg sim -r -t ${CPU} -x ${GRAPHS_PREFIX}.xg -F ${READS_PREFIX}_subset_1.fq -F ${READS_PREFIX}_subset_2.fq -T ${PROFILE_PREFIX}_haplotype.isoforms.results -n ${NREADS} -s ${SEED} -d 0.001 -p 216 -v 24 -a -E ${OUT_PREFIX}.txt > ${OUT_PREFIX}.gam"

# Extract simulated reads from alignments
/usr/bin/time -v bash -c "vg view -a -X ${OUT_PREFIX}.gam > ${OUT_PREFIX}.fq; wc -l ${OUT_PREFIX}.fq"

# De-interleave simulated reads (https://gist.github.com/nathanhaigh/3521724)
/usr/bin/time -v bash -c 'cat '"${OUT_PREFIX}"'.fq | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > '"${OUT_PREFIX}"'_1.fq) | cut -f 5-8 | tr "\t" "\n" > '"${OUT_PREFIX}"'_2.fq; wc -l '"${OUT_PREFIX}"'_1.fq; wc -l '"${OUT_PREFIX}"'_2.fq'

# Compress simulated reads and read info
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}_1.fq; gzip ${OUT_PREFIX}_2.fq; gzip ${OUT_PREFIX}.txt"

# Upload simulated reads and alignments
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/${TRANSCRIPTS_PREFIX}/${EXPRESSION_PREFIX}/vg/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "${OUT_PREFIX}.fq" --exclude "${OUT_PREFIX}.gam" --no-progress
