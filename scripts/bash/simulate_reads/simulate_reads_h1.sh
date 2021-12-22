set -e

# Set rng seed
SEED=7746153

# Set fragment mean and sd
MEAN=216
SD=24

# Set read (FASTQ) files
READ_1="ENCFF001REK.fastq.gz"
READ_2="ENCFF001REJ.fastq.gz"

# Set graph (XG) prefix
GRAPH_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set transcript info (txt) prefix
TRANSCRIPTS_PREFIX="1kg_NA12878_exons_gencode100_allpaths"

# Set expression profile (isoforms.results) prefix
PROFILE_PREFIX="1kg_NA12878_gencode100_ENCSR000AED_rep1_uni"

# Set output name prefix
OUT_PREFIX="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1"

# Set number of threads
CPU=1

# Extract haplotype-specific expression
/usr/bin/time -v bash -c "zcat ${TRANSCRIPTS_PREFIX}.txt.gz | grep NA12878 | grep _0_ | cut -f1 | grep -v Name > haplotype_transcripts.txt; echo transcript_id >> haplotype_transcripts.txt; wc -l haplotype_transcripts.txt; grep -F -f haplotype_transcripts.txt ${PROFILE_PREFIX}.isoforms.results > expression_haplotype.isoforms.results; wc -l expression_haplotype.isoforms.results; rm haplotype_transcripts.txt"

# Subset read pairs
/usr/bin/time -v bash -c "seqtk sample -s ${SEED} ${READ_1} 10000000 > subset_1.fq; seqtk sample -s ${SEED} ${READ_2} 10000000 > subset_2.fq; wc -l subset_1.fq; wc -l subset_2.fq"

# Simulate reads
/usr/bin/time -v bash -c "vg sim -r -t ${CPU} -x ${GRAPHS_PREFIX}.xg -F subset_1.fq -F subset_2.fq -T expression_haplotype.isoforms.results -n 25000000 -s ${SEED} -d 0.001 -p ${MEAN} -v ${SD} -a -E ${OUT_PREFIX}.txt > ${OUT_PREFIX}.gam; rm expression_haplotype.isoforms.results; rm subset*"

# Extract simulated reads from alignments
/usr/bin/time -v bash -c "vg view -a -X ${OUT_PREFIX}.gam > ${OUT_PREFIX}.fq; wc -l ${OUT_PREFIX}.fq"

# De-interleave simulated reads (https://gist.github.com/nathanhaigh/3521724)
/usr/bin/time -v bash -c 'cat '"${OUT_PREFIX}"'.fq | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > '"${OUT_PREFIX}"'_1.fq) | cut -f 5-8 | tr "\t" "\n" > '"${OUT_PREFIX}"'_2.fq; wc -l '"${OUT_PREFIX}"'_1.fq; wc -l '"${OUT_PREFIX}"'_2.fq; rm '"${OUT_PREFIX}"'.fq'

# Compress simulated reads and read info
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}_1.fq; gzip ${OUT_PREFIX}_2.fq; gzip ${OUT_PREFIX}.txt"
