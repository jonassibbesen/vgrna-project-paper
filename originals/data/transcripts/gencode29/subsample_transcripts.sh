grep -P "\ttranscript\t" gencode.v29.primary_assembly.annotation_renamed_full.gtf | cut -f9 | cut -d ';' -f2 | cut -d '"' -f2 | uniq | shuf | head -n $(echo "172449 * ${1} / 100" | bc | cut -d '.' -f1) > transcripts_subset${1}.txt

grep -F -f <(sed -e 's/$/\";/g' transcripts_subset${1}.txt) gencode.v29.primary_assembly.annotation_renamed_full.gtf > gencode.v29.primary_assembly.annotation_renamed_full_subset${1}.gtf

grep -F -f <(sed -e 's/$/;/g' transcripts_subset${1}.txt | sed -e 's/^/Parent=/g') gencode.v29.primary_assembly.annotation_renamed_full.gff3 > gencode.v29.primary_assembly.annotation_renamed_full_subset${1}.gff3
