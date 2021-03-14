## Construct spliced pangenome graph and pantranscriptome

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Construct graphs

`construct_graph.sh`: Construct spliced pangenome graph for a chromosome. Run this script on each chromosome including mitochondria and the scaffolds (`seq 1 22; echo "X"; echo "Y"; echo "MT"; echo "SCA"`).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

Important variables:

* `CHR`: Chromosome name

  

### 2. Join graphs

`join_graphs.sh`: Join the id space of the chromosome graphs. 

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1



### 3. Project transcripts

`project_transcripts.sh`: Project transcripts onto haplotypes and create pantranscriptome (haplotype-specifc trasncripts) for a chromosome. Run this script on each chromosome (`seq 1 22; echo "X"; echo "Y"`). Use `project_transcripts_gene.sh` instead if an exon-only splice graph is wanted.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

Important variables:

* `MODE`: Type of pantranscriptome to create (see script for more details)
* `CHR`: Chromosome name


### 4. Generate xg graph

`generate_xg.sh`: Combine chromosome graphs and create a xg graph for the whole genome

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

  

### 5. Merge pantranscriptomes

`merge.sh`: Merge chromosome pantranscriptomes (GBWT indexes, FASTA sequences and haplotype-specifc transcript info files).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1



### 6. Create r-index

`create_rindex.sh`: Create r-index of the merged GBWT index.

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1



### (7. Calculate pantranscriptome stats)

`calculate_hst_stats.sh`: Calculate population uniqueness statistics for haplotype-specific transcripts in the pantranscriptome. Use this population information file for 1000 Genomes Project data: https://ftp-trace.ncbi.nlm.nih.gov/1000genomes//ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2

