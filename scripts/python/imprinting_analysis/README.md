# Imprinting analyses

Scripts used for the analysis of imprinting in the vg rna/vg mpmap/rpvg paper.

## Reproducing the analysis

These scripts have the following dependencies:

* Executables `bcftools`, `tabix`, and `bgzip` discoverable on `$PATH`
* Python packages `pandas`, `numpy`, `matplotlib`, and `seaborn`

And the following data is required as input:

* VCF files used to construct the pantranscriptome (listed as in `example_input/vcf_list.txt`)
* List of focal imprinted genes (`example_input/focal_genes.txt` contains the imprinted genes used in the paper)
* List of samples to analyze (listed as in `example_input/sample_list.txt`)
* Expression output from `rpvg` (`expression_XXXXXX.txt`, not included in repo)
* Expression Gibbs sampling output from `rpvg` (`gibbs_XXXXXX.txt`, not included in repo)
* Gene annotation used to construct pantranscriptome (`annotation.gtf`, not included in repo)
* Allele to haplotype-specific transcript tables from `construct_allele_hst_table.py` (listed as in `example_input/hst_table_list.txt`)

The following pipeline will then reproduce the paper's analysis:

	mkdir sample_vcfs
	./preprocess_vcfs.py \
		annotation.gtf \
		example_input/focal_genes.txt \
		example_input/vcf_list.txt \
		example_input/sample_list.txt \
		sample_vcfs > genotypes.tsv
	
	for f in `ls sample_vcfs | grep -v ".tbi$"`; do realpath sample_vcfs/$f >> sample_vcf_list.txt; done
	
	mkdir imprinting_output
	./imprinting_analysis.py \
		annotation.gtf \
		example_input/focal_genes.txt \
		NA12878:expression_NA12878.txt,NA12891:expression_NA12891.txt \
		NA12878:gibbs_NA12878.txt,NA12891:gibbs_NA12891.txt \
		example_input/hst_table_list.txt \
		sample_vcf_list.txt \
		genotypes.tsv \
		imprinting_output
	
	mkdir plots
	./make_plots.py \
		imprinting_output \
		example_input/focal_genes.txt \
		plots
