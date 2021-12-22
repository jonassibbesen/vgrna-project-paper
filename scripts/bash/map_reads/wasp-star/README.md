## Map RNA-seq reads using STAR-WASP pipeline

The script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 

### 1. Map reads

`map_reads.sh`: Map paired-end RNA-seq reads to spliced reference genome and filter allele-biased reads with WASP.

* Docker image: quay.io/jsibbesen/star-s3script:star-2.7.3a-s1


