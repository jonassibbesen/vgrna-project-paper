set -e

# Set mapping quality threshold
MAPQ=30

# Set output name prefix
OUT_PREFIX="ENCSR706ANY_mq${MAPQ}"

# Set number of threads
CPU=1

# Download alignments
aws s3 cp s3://encode-public/2019/11/20/79bacc96-a43f-4130-bfc3-1d05550545d2/ENCFF247TLH.bam . --no-sign-request --no-progress
aws s3 cp s3://encode-public/2019/11/20/fa2356af-4075-4fe8-bebc-553b71444d0f/ENCFF431IOE.bam . --no-sign-request --no-progress
aws s3 cp s3://encode-public/2019/11/20/b2de9f99-6c82-46bb-bf70-9cd535fc52a9/ENCFF520MMC.bam . --no-sign-request --no-progress
aws s3 cp s3://encode-public/2019/11/20/0363ca4c-0790-4c38-8829-f87341ec6129/ENCFF626GWM.bam . --no-sign-request --no-progress

# Merge alignments
/usr/bin/time -v bash -c "samtools merge -O BAM merged.bam ENCFF247TLH.bam ENCFF431IOE.bam ENCFF520MMC.bam ENCFF626GWM.bam"

# Update alignment header
/usr/bin/time -v bash -c "samtools view -H  merged.bam | sed -e 's/SN:chr/SN:/g' | sed -e 's/SN:M/SN:MT/g' | sed -e 's/_random//g' | sed -E 's/SN:.*_(.*)/SN:\1/g' | sed -E 's/v([1-2]{1})/.\1/g' > new_header.sam; samtools reheader new_header.sam merged.bam > merged_header.bam; rm new_header.sam"

# Filter alignments
/usr/bin/time -v bash -c "samtools view -F 256 -q ${MAPQ} -O BAM merged_header.bam > merged_header_filter.bam"

# Sort alignments
/usr/bin/time -v bash -c "samtools sort -O BAM merged_header_filter.bam > ${OUT_PREFIX}.bam; samtools index ${OUT_PREFIX}.bam; samtools flagstat ${OUT_PREFIX}.bam; rm merged*"

# Find alignment exons
/usr/bin/time -v bash -c "bedtools bamtobed -splitD -i ${OUT_PREFIX}.bam | cut -f1-3 > exons.bed; bedtools sort -i exons.bed > exons_sort.bed; bedtools merge -i exons_sort.bed > ${OUT_PREFIX}.bed; wc -l ${OUT_PREFIX}.bed; rm exons.bed"
