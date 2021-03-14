## Map RNA-seq reads using HISAT2

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Prepare variants

`prepare_variants.sh`: Construct variant and haplotype list for a chromosome. Run this script on each chromosome (`seq 1 22; echo "X"; echo "Y"`).

* Docker image: quay.io/jsibbesen/hisat2-s3script:hisat2-2.2.1-s2

Important variables:

* `CHR`: Chromosome name



### 2. Generate index

`generate_index.sh`: Generate index of spliced pangenome graph.

* Docker image: quay.io/jsibbesen/hisat2-s3script:hisat2-2.2.1-s2



### 3. Map reads

`map_reads.sh`: Map paired-end RNA-seq reads to spliced pangenome graph.

* Docker image: quay.io/jsibbesen/hisat2-s3script:hisat2-2.2.1-s2



### (4. Inject alignments)

`inject_alignments.sh`: Optionally inject alignments to graph coordinates (needed for graph-based mapping evaluation).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

