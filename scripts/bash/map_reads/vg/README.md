## Map RNA-seq reads using vg map or vg mpmap

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Generate distance index

`generate_distance.sh`: Generate distance index of spliced pangenome graph.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1



### 2. Generate GCSA index

`generate_index.sh`: Generate GCSA index of spliced pangenome graph.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1



### 3. Map reads

`map_reads.sh`: Map paired-end RNA-seq reads to spliced pangenome graph.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-2cea1e25-s2

Important variables:

* `MODE`: Alignment algorithm to use (see script for more details)



### (4. Surject alignments)

`surject_alignments.sh`: Optionally surject alignments to reference coordinates (needed for reference-based mapping evaluation).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-2cea1e25-s2

