grep -P "^${1}\t" /public/groups/vg/jsibbese/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full.gtf > gencode_${1}.gtf

paste <(grep -P "\texon\t" gencode_${1}.gtf | cut -f1,4-5) <(grep -P "\texon\t" gencode_${1}.gtf | cut -f9 | cut -d ';' -f2 | cut -d '"' -f2) | sort -k2 -k3 -n | bgzip > gencode_${1}.txt.gz
tabix -s1 -b2 -e3 gencode_${1}.txt.gz

/public/groups/vg/jsibbese/vgrna/tools/bcftools-1.11/bcftools annotate -a gencode_${1}.txt.gz -h ../header.hdr -c CHROM,FROM,TO,TRANSCIPTS -l TRANSCIPTS:unique 1kg_NA12878_exons_${1}.vcf.gz | gzip > 1kg_NA12878_exons_anno_${1}.vcf.gz

zcat 1kg_NA12878_exons_gencode100_allpaths_${1}.txt.gz | grep -E "Haplotypes|NA12878" | gzip -c > 1kg_NA12878_exons_gencode100_allpaths_${1}_noref.txt.gz

python /public/groups/vg/jsibbese/vgrna/code/vgrna-project-scripts/python/construct_allele_hst_table.py 1kg_NA12878_exons_anno_${1}.vcf.gz 1kg_NA12878_exons_gencode100_allpaths_${1}_noref.txt.gz hst_variants_${1}.txt
gzip hst_variants_${1}.txt
