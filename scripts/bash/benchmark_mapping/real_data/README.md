## Benchmark RNA-seq mapping using real data

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### (1. Calculate Iso-Seq exon coverage statistics)

Optionally calculate exon coverage for Iso-Seq alignments from the same cell-line (see [ENCSR706ANY](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/benchmark_mapping/real_data/ENCSR706ANY) folder).



### 1. Calculate overlap and coverage statistics

`calculate_bam_stats.sh`: Calculate exon overlap and coverage using an annotation and optionally exon regions defined by Iso-Seq alignments.

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-2220bb08-s7



### 2. Plot benchmark

All real data mapping benchmark plots presented in the paper were created using the following R scripts from [here](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/scripts/R):

* plot_mapping_stats.R
* plot_coverage_correlation.R
* plot_mapping_compute.R
* plot_mapping_memory.R

Please note that some of the scripts assume that the data is structured in a specific set of folders. They might therefore require a bit of reworking before they can be used in other environments. If there is an interest for it we can look into creating more user-friendly versions of the scripts. 
