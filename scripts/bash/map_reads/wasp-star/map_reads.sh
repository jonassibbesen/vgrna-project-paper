set -e

# Set read (FASTQ) files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_2.fq.gz"
# Set index directory
INDEX_DIR="gencode100_index"
# Set output name prefix
OUT_PREFIX="star_wasp_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

WASP_DIR="~/GitHub/WASP/"

echo "Initial mapping with STAR"
/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir ${INDEX_DIR} --readFilesCommand zcat --readFilesIn ${READ_1} ${READ_2} --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_; mv ${OUT_PREFIX}_Aligned.out.bam ${OUT_PREFIX}.bam"

echo "Flagging mapped reads with WASP"
/usr/bin/time -v bash -c "python3 ${WASP_DIR}/mapping/find_intersecting_snps.py --is_paired_end --output_dir . --snp_index snp_index.h5 --snp_tab snp_tab.h5 --haplotype haplotypes.h5 --sample ${POPULATION}_samples.txt ${OUT_PREFIX}.bam"

echo "Re-mapping with STAR"
/usr/bin/time -v bash -c "STAR --runThreadN ${CPU} --genomeDir ${INDEX_DIR} --readFilesCommand zcat --readFilesIn ${OUT_PREFIX}.remap.fq1.gz ${OUT_PREFIX}.remap.fq2.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${OUT_PREFIX}_remap_; mv ${OUT_PREFIX}_remap_Aligned.out.bam ${OUT_PREFIX}_remap.bam"

echo "Filtering reads with WASP"
/usr/bin/time -v bash -c "python3 ${WASP_DIR}/mapping/filter_remapped_reads.py ${OUT_PREFIX}.to.remap.bam ${OUT_PREFIX}_remap.bam ${OUT_PREFIX}_remap_filtered.bam"


echo "Merging filtered and unfiltered reads"
/usr/bin/time -v bash -c "samtools merge ${OUT_PREFIX}_filtered.bam ${OUT_PREFIX}_remap_filtered.bam ${OUT_PREFIX}.keep.bam"

echo "Sorting and indexing"
/usr/bin/time -v bash -c "samtools sort ${OUT_PREFIX}_filtered.bam -o ${OUT_PREFIX}_filtered_sorted.bam; mv ${OUT_PREFIX}_filtered_sorted.bam ${OUT_PREFIX}_filtered.bam; samtools index ${OUT_PREFIX}_filtered.bam"

echo "Removing duplicates with WASP"
/usr/bin/time -v bash -c "python3 ${WASP_DIR}/mapping/rmdup_pe.py ${OUT_PREFIX}_filtered.bam ${OUT_PREFIX}_final.bam"


