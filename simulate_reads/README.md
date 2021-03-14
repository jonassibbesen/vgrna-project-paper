## Simulate RNA-seq reads

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Construct transcripts

`construct_graph.sh`: Construct graph and haplotype-specific transcripts for a single sample using pipeline described under [contruct_pantranscriptome](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/contruct_pantranscriptome).



### 2. Infer expression

`infer_expression.sh`: Infer expression profile of haplotype-specific transcripts.

* Docker image: quay.io/jsibbesen/rsem-s3script:rsem-1.3.1-s1

Important variables:

* `SEED`: Seed for the random number generator



### (2+. Uniform expression)

`uniform_expression.sh`: Optionally convert expression profile to uniform expression across all transcripts.

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2



### 3. Simulate reads

`simulate_reads_h*.sh`: Simulate paired-end RNA-seq reads from haplotype-specific transcripts on each haplotype separately (`h1` and `h2`).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1

Important variables:

* `SEED`: Seed for the random number generator
* `MEAN`: Mean of fragment length distribution
* `SD`: Standard deviation of fragment length distribution



### (4. Surject transcripts)

`surject_transcripts.sh`: Optionally surject haplotype-specific transcripts to reference coordinates (needed for reference-based mapping evaluation).

* Docker image: quay.io/jsibbesen/vgdev-s3script:vgdev-c4bbd63b-s1