## Infer haplotype-specific transcript expression using rpvg

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Infer expression

`infer_expression.sh`: Infer expression of haplotype-specific transcripts in pantranscriptome.

* Docker image: quay.io/jsibbesen/rpvg-s3script:rpvg-a7a79697-s9

Important variables:

* `MODE`: Inference mode to use (see script for more details)
* `SEED`: Seed for the random number generator
* `MEAN`: Mean of fragment length distribution (only used in single-path mode)
* `SD`: Standard deviation of fragment length distribution (only used in single-path mode)
