for i in $(seq 1 22; echo "X"; echo "Y"); do zcat ${i}/1kg_all_exons_${i}_vep.txt.gz | cut -f1-16,21,25,31,38-43,48,56,67-69,71-75 | gzip > ${i}/1kg_all_exons_${i}_vep_subset.txt.gz; done
