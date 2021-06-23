set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes"
TRANSCRIPTS_PREFIX="gencode.v29.primary_assembly.annotation_renamed_full"
BASE_PREFIX="1kg_nonCEU_af001_gencode100"
GRAPHS_PREFIX="1kg_nonCEU_af001_imgt_hla_noN_main_gencode100"
OUT_PREFIX="${GRAPHS_PREFIX}_6"

# Download genome
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa genome.fa --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/genomes/GRCh38/${GENOME_PREFIX}.fa.fai genome.fa.fai --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS_PREFIX}.gtf transcripts.gtf --no-progress

# Download graphs
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${BASE_PREFIX}/${BASE_PREFIX}.xg base.xg --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/${GRAPHS_PREFIX}_6.pg 6.pg --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${BASE_PREFIX}/${BASE_PREFIX}_index.gcsa index.gcsa --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${BASE_PREFIX}/${BASE_PREFIX}_index.gcsa.lcp index.gcsa.lcp --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${BASE_PREFIX}/${BASE_PREFIX}_index.dist index.dist --no-progress

# Download HLA alleles
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/haplotypes/imgt_hla/A/A_cds_alleles_full_pad1k.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/haplotypes/imgt_hla/B/B_cds_alleles_full_pad1k.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/haplotypes/imgt_hla/C/C_cds_alleles_full_pad1k.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/haplotypes/imgt_hla/DQB1/DQB1_cds_alleles_full_pad1k.fa . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/haplotypes/imgt_hla/DRB1/DRB1_cds_alleles_full_pad1k.fa . --no-progress

# Filter HLA CDS null alleles
/usr/bin/time -v bash -c "wc -l *_cds_alleles_full_pad1k.fa; grep '>' *_cds_alleles_full_pad1k.fa | wc -l; grep -h '>' *_cds_alleles_full_pad1k.fa | grep -v 'N$' | sed -e 's/>//g' > no_null_alleles.txt; wc -l no_null_alleles.txt; seqtk subseq <(cat *_cds_alleles_full_pad1k.fa) no_null_alleles.txt > cds_alleles_full_pad1k_non.fa; wc -l cds_alleles_full_pad1k_non.fa; grep '>' cds_alleles_full_pad1k_non.fa | wc -l"

# Map HLA CDS alleles
/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -l long -F GAM -x base.xg -g index.gcsa -d index.dist -f cds_alleles_full_pad1k_non.fa > cds.gam"

# Surject CDS alignments
/usr/bin/time -v bash -c "vg surject -t ${CPU} -S -b -p 6 -x base.xg cds.gam > ${OUT_PREFIX}_cds.bam"

# Calculate CDS alignment statistics
/usr/bin/time -v bash -c "samtools flagstat ${OUT_PREFIX}_cds.bam"

# Find contig transcripts
/usr/bin/time -v bash -c "grep -P '^6\t' transcripts.gtf > 6.gtf"

# Convert CDS alignments to haplotypes
/usr/bin/time -v bash -c "convert_cds_alignments_to_haplotypes ${OUT_PREFIX}_cds.bam genome.fa 6.gtf 2 | sed -e 's/>/>hla_/g' > ${OUT_PREFIX}_haps.fa"

# Map HLA haplotypes
/usr/bin/time -v bash -c "grep '>' ${OUT_PREFIX}_haps.fa | wc -l; vg mpmap -t ${CPU} -l long -F GAM -x base.xg -g index.gcsa -d index.dist -f ${OUT_PREFIX}_haps.fa > haps.gam"

# Calculate haplotype alignment statistics
/usr/bin/time -v bash -c "vg stats -a haps.gam"

# Calculate graph statistics
/usr/bin/time -v bash -c "vg stats -z -l -r 6.pg"

# Augment graph and haplotypes alignments
/usr/bin/time -v bash -c "vg augment -p -t ${CPU} -i 6.pg haps.gam | vg convert -t ${CPU} -p - > ${OUT_PREFIX}.pg"

# Calculate updated graph statistics
/usr/bin/time -v bash -c "vg stats -z -l -r ${OUT_PREFIX}.pg"

# Upload updated graph and haplotypes
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPHS_PREFIX}/6/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
