set -e

OUT_PREFIX="STAR_diploid_real_r1_ENCSR000AED_rep1"
SCRIPTS_DIR="~/GitHub/vgrna-project-paper/scripts/"
CPU=1
INDEX_DIR_SUFFIX="_indexes"
MATERNAL_INDEX_DIR="maternal${INDEX_DIR_SUFFIX}"
PATERNAL_INDEX_DIR="paternal${INDEX_DIR_SUFFIX}"
sudo pip3 install CrossMap

echo "Making joint reverse chain"

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap
chmod +x chainSwap
chainSwap ${MATERNAL_INDEX_DIR}/maternal.chain maternal.reverse.chain
chainSwap ${PATERNAL_INDEX_DIR}/paternal.chain paternal.reverse.chain

cat *.reverse.chain > joint.reverse.chain

echo "Mapping reads"

# map to each haplotype independently
for PARENTAL in maternal paternal; do

    echo "Mapping reads to $PARENTAL indexes"
    
    /usr/bin/time -v ./STAR-2.7.9a/bin/Linux_x86_64_static/STAR --runThreadN ${CPU} --genomeDir ${PARENTAL}${INDEX_DIR_SUFFIX} --readFilesCommand zcat --readFilesIn reads_1.fq.gz reads_2.fq.gz --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix ${PARENTAL}_${OUT_PREFIX}_
    mv ${PARENTAL}_${OUT_PREFIX}_Aligned.out.bam ${PARENTAL}_${OUT_PREFIX}.bam
    
    echo "Adding read number to read name"
    
    # add the read 2 back onto the read name
    /usr/bin/time -v samtools merge -O BAM ${PARENTAL}_${OUT_PREFIX}.numbered.bam <(samtools view -h -f 128 ${PARENTAL}_${OUT_PREFIX}.bam | awk -v OFS='\t' '{ if (substr($1, 1, 1) != "@") gsub("$", "_2", $1) } 1' | samtools view -O BAM ) <(samtools view -h -F 128 ${PARENTAL}_${OUT_PREFIX}.bam | awk -v OFS='\t' '{ if (substr($1, 1, 1) != "@") gsub("$", "_1", $1)} 1' | samtools view -O BAM )
    rm ${PARENTAL}_${OUT_PREFIX}.bam
    
    echo "Filtering to primary alignments and sorting by name"
    
    # filter out secondary alignments and sort by name to make deterministic order for next step
    /usr/bin/time -v samtools sort -n -@ ${CPU} -O BAM <(samtools view -b -F 256 ${PARENTAL}_${OUT_PREFIX}.numbered.bam) > ${PARENTAL}_${OUT_PREFIX}.namesorted.bam
    rm ${PARENTAL}_${OUT_PREFIX}.numbered.bam
    
    echo "Removing read number from read names"
    
    /usr/bin/time -v samtools view -h ${PARENTAL}_${OUT_PREFIX}.namesorted.bam | awk -v OFS='\t' '{gsub("_[12]$", "", $1)} 1' | samtools view -@ ${CPU} -b > ${PARENTAL}_${OUT_PREFIX}.restorename.bam
    rm ${PARENTAL}_${OUT_PREFIX}.namesorted.bam
done
    
echo "Choosing haplotypes for reads"

# choose a haplotype for each read
/usr/bin/time -v ${SCRIPTS_DIR}/python/choose_mapping.py <(samtools view -h maternal_${OUT_PREFIX}.restorename.bam) <(samtools view -h paternal_${OUT_PREFIX}.restorename.bam) | samtools view -b -@ ${CPU} - > joint_${OUT_PREFIX}.bam

/usr/bin/time -v samtools sort -@ ${CPU} joint_${OUT_PREFIX}.bam > joint_${OUT_PREFIX}.sorted.bam
rm joint_${OUT_PREFIX}.bam
samtools index joint_${OUT_PREFIX}.sorted.bam

echo "Lifting to reference"

# lift back onto the reference
/usr/bin/time -v CrossMap.py bam --mean 300 --stdev 75 joint.reverse.chain joint_${OUT_PREFIX}.sorted.bam ${OUT_PREFIX} # adds .sorted.bam on its own
rm joint_${OUT_PREFIX}.sorted.bam
samtools view -@ ${CPU} -h ${OUT_PREFIX}.sorted.bam | ${SCRIPTS_DIR}/python/fix_cross_map_flags.py | samtools view -b -@ ${CPU} - > ${OUT_PREFIX}.bam
rm ${OUT_PREFIX}.sorted.bam

