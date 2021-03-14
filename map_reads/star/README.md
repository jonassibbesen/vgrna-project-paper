## Map RNA-seq reads using STAR

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Generate index

`generate_index.sh`: Generate index of spliced reference genome.

* Docker image: quay.io/jsibbesen/star-s3script:star-2.7.3a-s1



### 2. Map reads

`map_reads.sh`: Map paired-end RNA-seq reads to spliced reference genome.

* Docker image: quay.io/jsibbesen/star-s3script:star-2.7.3a-s1



### (3. Inject alignments)

`inject_alignments.sh`: Optionally inject alignments to graph coordinates (needed for graph-based mapping evaluation).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

