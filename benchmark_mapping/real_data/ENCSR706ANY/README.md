## Calculate Iso-Seq alignment (ENCSR706ANY) exon coverage

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Prepare alignments

`prepare_alignments.sh`: Merge, rename contigs and filter ENCSR706ANY alignments. Create exon regions defined by the alignments.

* Docker image: jsibbesen/base-s3script:18.04-s1

Important variables:

* `MAPQ`: Mapping quality thresholds



### 2. Calculate exon coverage

`calculate_exon_coverage.sh`: Calculate exon coverage in both regions defined by the alignments and an annotation.

* Docker image: quay.io/jsibbesen/vgrna-s3script:vgrna-71442ea4-s2
