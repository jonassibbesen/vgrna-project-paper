# Haplotype-aware pantranscriptome analyses using spliced pangenome graphs

This repository contains the scripts that were used to generate the results presented in *Haplotype-aware pantranscriptome analyses using spliced pangenome graphs*.  

For more up-to-date information on how to run the different methods please go to github pages of the [vg toolkit](https://github.com/vgteam/vg) and [rpvg](https://github.com/jonassibbesen/rpvg). The spliced pangenome graphs and pantranscriptomes (haplotype-specific transcripts) presented in the paper are avaliable to download in the [Data](#Data) section for use in other projects.

There are two types of script in this repository: 

1. In the [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) folder we include the raw scripts and the log files. These files are not particularly user-friendly as they include a lot of hard-coded paths, however, we have included them here for transparency and reproducibility. By looking at the scripts and log files you can see exactly how each method was run in the paper. Most of the log files will include a short header which specifies the Docker image that was used. The Docker containers are available at this repository: [s3script-dockerfiles](https://github.com/jonassibbesen/s3script-dockerfiles). For the log files without this header it should be clear from the script itself what version was used. All plotting was done using the R scripts in the [vgrna-project-scripts](https://github.com/jonassibbesen/vgrna-project-scripts) repository. 

2. The remaining folders include versions of the scripts in [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) that have been updated to make them easier to be used by others. These changes mostly includes removing hard-coded paths and adding variables for the input and output files that can be set in the header by the user. These updated scripts should in theory be able to produce the same results as presented in the paper, however, we can not guarantee that (the scripts in [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) can). Each folder contains a short readme that explains each script, their order and which Docker images to use. 

## Data

Here you can find links to the data used in the paper. This includes both the raw data and data constructed as part of the analyses in the paper. The constructed data included here is either the data that are not guaranteed to be reproducible (subsampled transcript annotation and simulated reads) or that are deemed potentially useful in other projects (graphs, pantranscriptomes and indexes).

### Genome

The GRCh38 (primary assembly) reference genome:

* ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/

### Transcripts

The GENCODE v29 (primary assembly) transcript annotation:

* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz 

The subsampled (80%) GENCODE v29 transcript annotation:

* http://cgl.gi.ucsc.edu/data/vgrna/transcript_annotation/

### Variants

The 1000 Genomes Project variants and haplotypes lifted to GRCh38:

* http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ 

### Reads

The simulated RNA-seq reads:

* http://cgl.gi.ucsc.edu/data/vgrna/simulated_data/

The real RNA-seq reads:

* https://www.ncbi.nlm.nih.gov/sra/?term=SRR1153470
* https://www.encodeproject.org/experiments/ENCSR000AED/ (replicate 1)
* https://github.com/nanopore-wgs-consortium/CHM13/ (replicate 1)

The Iso-Seq alignments:

* https://www.encodeproject.org/experiments/ENCSR706ANY/ (all replicates)

### Graphs, pantranscriptomes and indexes

The spliced pangenome graphs, pantranscriptomes and indexes:

* http://cgl.gi.ucsc.edu/data/vgrna/pantranscriptomes/
