zcat gencode.v29.primary_assembly.annotation.gff3.gz | sed -e 's/^chrM/chrMT/g' | grep -v "mRNA_start_NF" | grep -v "mRNA_end_NF" > gencode.v29.primary_assembly.annotation_renamed_full.gff3
sed -i -e 's/^chr//g' gencode.v29.primary_assembly.annotation_renamed_full.gff3
zcat gencode.v29.primary_assembly.annotation.gtf.gz | sed -e 's/^chrM/chrMT/g' | grep -v "mRNA_start_NF" | grep -v "mRNA_end_NF" > gencode.v29.primary_assembly.annotation_renamed_full.gtf
sed -i -e 's/^chr//g' gencode.v29.primary_assembly.annotation_renamed_full.gtf
