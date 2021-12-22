eset -e

echo "Setting up"

# get dependencies
sudo apt update
sudo apt install -qq -y python3
python3 -c 'print("successfully installed python")'
sudo apt install -qq -y time awscli git-all python3 python3-pip samtools
sudo pip install numpy pysam tables

# get STAR
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz
tar -xzf 2.7.9a.tar.gz
sudo cp STAR-2.7.9a/bin/Linux_x86_64_static/STAR /usr/local/bin

# get WASP
git clone https://github.com/bmvdgeijn/WASP.git

echo "Downloading files"

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/star/indexes/${REF}/ . --recursive --no-progress
aws s3 cp s3://vg-k8s/users/jeizenga/wasp/h5${VARS}/ . --recursive --no-progress
aws s3 cp s3://vg-k8s/users/jeizenga/wasp/samples/${POPULATION}_samples.txt . --no-progress

# get reads
aws s3 cp ${READS_1} reads_1.fq.gz --no-progress
aws s3 cp ${READS_2} reads_2.fq.gz --no-progress

# Set file name prefixes
OUT_PREFIX="${MAPPER}_${REF}_${TYPE}_${REAL}"

echo "Initial mapping with STAR"
/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir . --readFilesCommand zcat --readFilesIn reads_1.fq.gz reads_2.fq.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_; mv ${OUT_PREFIX}_Aligned.out.bam ${OUT_PREFIX}.bam"

echo "Flagging mapped reads with WASP"
/usr/bin/time -v bash -c "python3 ./WASP/mapping/find_intersecting_snps.py --is_paired_end --output_dir . --snp_index snp_index.h5 --snp_tab snp_tab.h5 --haplotype haplotypes.h5 --sample ${POPULATION}_samples.txt ${OUT_PREFIX}.bam"

echo "Re-mapping with STAR"
/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir . --readFilesCommand zcat --readFilesIn ${OUT_PREFIX}.remap.fq1.gz ${OUT_PREFIX}.remap.fq2.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_remap_; mv ${OUT_PREFIX}_remap_Aligned.out.bam ${OUT_PREFIX}_remap.bam"

echo "Filtering reads with WASP"
/usr/bin/time -v bash -c "python3 ./WASP/mapping/filter_remapped_reads.py ${OUT_PREFIX}.to.remap.bam ${OUT_PREFIX}_remap.bam ${OUT_PREFIX}_remap_filtered.bam"


echo "Merging filtered and unfiltered reads"
/usr/bin/time -v bash -c "samtools merge ${OUT_PREFIX}_filtered.bam ${OUT_PREFIX}_remap_filtered.bam ${OUT_PREFIX}.keep.bam"

echo "Sorting and indexing"
/usr/bin/time -v bash -c "samtools sort ${OUT_PREFIX}_filtered.bam -o ${OUT_PREFIX}_filtered_sorted.bam; mv ${OUT_PREFIX}_filtered_sorted.bam ${OUT_PREFIX}_filtered.bam; samtools index ${OUT_PREFIX}_filtered.bam"

echo "Removing duplicates with WASP"
/usr/bin/time -v bash -c "python3 ./WASP/mapping/rmdup_pe.py ${OUT_PREFIX}_filtered.bam ${OUT_PREFIX}_final.bam"

# Upload read alignments
aws s3 cp ${OUT_PREFIX}_final.bam s3://vg-k8s/users/jeizenga/wasp/${OUT_PREFIX}_`date +"%N"`.bam --no-progress

