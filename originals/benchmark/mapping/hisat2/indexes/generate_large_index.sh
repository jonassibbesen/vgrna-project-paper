set -e

# Set file name prefixes 
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"
OUT_PREFIX="${REF}_index"

# Download genome
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/${GENOME_PREFIX}.fa.fai . --no-progress

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/${TRANSCRIPTS} . --no-progress

# Construct exon list
/usr/bin/time -v bash -c "hisat2_extract_exons.py ${TRANSCRIPTS} > exons.txt; wc -l exons.txt"

# Construct splice-site list
/usr/bin/time -v bash -c "hisat2_extract_splice_sites.py ${TRANSCRIPTS} > splice_sites.txt; wc -l splice_sites.txt"

if [ "${VARIANTS}" = "" ]; then

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} --large-index --exon exons.txt --ss splice_sites.txt ${GENOME_PREFIX}.fa ${OUT_PREFIX}"

else

	# Download variants
	aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/variants/${VARIANTS}/ . --recursive --exclude "*" --include "*.snp" --include "*.haplotype" --no-progress

	# Combine variants and haplotype lists
	/usr/bin/time -v bash -c "cat */*.snp > variants.snp; cat */*.haplotype > variants.haplotype; wc -l variants.snp; wc -l variants.haplotype"

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} --large-index --snp variants.snp --haplotype variants.haplotype --exon exons.txt --ss splice_sites.txt ${GENOME_PREFIX}.fa ${OUT_PREFIX}"
fi

# Upload HISAT2 index
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/hisat2/indexes/${REF}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
