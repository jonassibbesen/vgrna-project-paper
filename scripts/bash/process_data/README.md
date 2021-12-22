## Processing the raw data

Each script includes a header with input and output file variables. These variables should be set before running the scripts. The number of threads can be set using `CPU`. 

Note that the Docker image suggested is not necessarily the same version of the container that was used in the paper. For information on which exact version was used see log files in the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder. 



### 1. Process genome

`process_genome.sh`: Split genome into two files containing chromosomes and scaffolds, respectively.

* Docker image: jsibbesen/base-s3script:18.04-s1

  

### 2. Process transcripts

`process_transcripts.sh`: Filter for full length transcripts, rename contigs and create exon region BED file. 

* Docker image: jsibbesen/base-s3script:18.04-s1



### 3. Process variants

`process_variants.sh`: Subset, filter and normalise variants on a chromosome. Run this script on each chromosome (`seq 1 22; echo "X"; echo "Y"`).

* Docker image: jsibbesen/base-s3script:18.04-s1

Important variables:

* `CHR`: Chromosome name
* `MAF`: Minimum allele frequency in intergenic regions
* `SAMPLES`: File containing the samples to subset



