## Benchmark RNA-seq mapping using simulated data

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Calculate reference-based statistics

`calculate_bam_stats.sh`: Calculate overlap between predicted alignments and simulated alignments on the reference-genome.

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2

Important variables:

* `MAPPER`: Mapper name



### 1. Calculate graph-based statistics

`calculate_gamp_stats.sh`: Calculate distance between predicted alignments and simulated alignments in the spliced pangenome graph.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

Important variables:

* `MAPPER`: Mapper name



### 1. Calculate bias statistics

`calculate_bam_bias.sh`: Calculate read alignment coverage across variants for each haplotype. Run this script on each chromosome (`seq 1 22; echo "X"; echo "Y"`).

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2

Important variables:

* `CHR`: Chromosome name



### 2. Plot benchmark

All simulated data mapping benchmark plots presented in the paper were created using the following R scripts from the [vgrna-project-scripts](https://github.com/jonassibbesen/vgrna-project-scripts) repository:

* plot_vg_sim_overlap.R
* plot_vg_sim_distance.R
* plot_mapping_bias.R

Please note that some of the scripts assume that the data is structured in a specific set of folders. They might therefore require a bit of reworking before they can be used in other environments. If there is an interest for it we can look into creating more user-friendly versions of the scripts. 
