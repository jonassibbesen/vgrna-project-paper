Rscript /public/groups/vg/jsibbese/vgrna/code/vgrna-project-scripts/R/impute_hla_cds_alleles.R ${1} ../database/IMGTHLA/ ${1}_cds_alleles_full.txt

python /public/groups/vg/jsibbese/vgrna/code/vgrna-project-scripts/python/reference_pad_cds_alleles.py ${1}_cds_alleles_full.txt /public/groups/vg/jsibbese/vgrna/data/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.fa /public/groups/vg/jsibbese/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation_renamed_full.gtf HLA-${1} 1000 ${1}_cds_alleles_full_pad1k.fa
