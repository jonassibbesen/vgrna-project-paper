## Benchmark haplotype-specific transcript expression inference

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Compare transcripts

`compare_hst.sh`: Compare haplotype-specific transcript sequences between two pantranscriptomes (e.g. compare a pantranscriptome to sequences used in simulation).

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2



### 2. Plot benchmark

All expression benchmark plots presented in the paper were created using the following R scripts from the [vgrna-project-scripts](https://github.com/jonassibbesen/vgrna-project-scripts) repository:

* parse_expression_data.R
* plot_expression_data.R
* plot_haploid_expression_data.R

Please note that some of the scripts assume that the data is structured in a specific set of folders. They might therefore require a bit of reworking before they can be used in other environments. If there is an interest for it we can look into creating more user-friendly versions of the scripts. 
