bcftools view -g het 1kg_NA12878_exons_${1}.vcf.gz | grep -v "^#" | cut -f1,2,4,5,10 > 1kg_NA12878_exons_${1}_het.txt
