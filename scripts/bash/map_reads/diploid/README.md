## Map RNA-seq reads using diploid reference STAR pipeline

The script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 

### 1. Generate indexes

`generate_index.sh`: Create maternal and paternal references using `vcf2diploid` and then index them for mapping by `STAR`.

### 2. Map reads

`map_reads.sh`: Ma paired-end RNA-seq reads to the diploid reference and choose a primary mapping for each read.

* Docker image: quay.io/jsibbesen/star-s3script


