python /public/groups/vg/jsibbese/vgrna/code/vgrna-project-scripts/python/calc_allele_rpvg_expression.py /public/groups/vg/jsibbese/vgrna/analyses/pantranscriptomes/1kg_all_af001_gencode100_unidi/${1}/1kg_all_exons_anno_${1}.vcf.gz /public/groups/vg/jsibbese/vgrna/analyses/pantranscriptomes/1kg_all_af001_gencode100_unidi/${1}/1kg_all_af001_gencode100_unidi_${1}.txt.gz /public/groups/vg/jsibbese/vgrna/data/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.fa ../../inference_r2/${2}/1kg_all_af001_gencode100_unidi/rpvg_mpmap_r2_${2}_1kg_all_af001_gencode100_unidi_joint.txt.gz ../../inference_r2/${2}/1kg_all_af001_gencode100_unidi/rpvg_mpmap_r2_${2}_1kg_all_af001_gencode100_unidi.txt.gz 200 allele_exp_rpvg_mpmap_r2_${2}_1kg_all_af001_gencode100_unidi_${1}.txt

gzip allele_exp_rpvg_mpmap_r2_${2}_1kg_all_af001_gencode100_unidi_${1}.txt
