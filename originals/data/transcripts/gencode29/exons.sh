grep -v "^#" gencode.v29.primary_assembly.annotation_renamed_full.gff3 | grep -P "\texon\t" | cut -f1,4,5 | awk -v OFS="\t" '{$2 -= 1}{print}' > gencode.v29.primary_assembly.annotation_renamed_full_exons.txt

/public/groups/vg/jsibbese/vgrna/tools/bedtools2/bin/bedtools sort -i gencode.v29.primary_assembly.annotation_renamed_full_exons.txt > gencode.v29.primary_assembly.annotation_renamed_full_exons_sort.txt
/public/groups/vg/jsibbese/vgrna/tools/bedtools2/bin/bedtools merge -i gencode.v29.primary_assembly.annotation_renamed_full_exons_sort.txt > gencode.v29.primary_assembly.annotation_renamed_full_exons.bed

rm gencode.v29.primary_assembly.annotation_renamed_full_exons*.txt
