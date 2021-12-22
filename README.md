# Haplotype-aware pantranscriptome analyses using spliced pangenome graphs

This repository contains the scripts that were used to generate the results presented in *Haplotype-aware pantranscriptome analyses using spliced pangenome graphs, [bioRxiv](https://doi.org/10.1101/2021.03.26.437240) (2021)*.  

For more up-to-date information on how to run the different methods please go to the github page of the [vg toolkit](https://github.com/vgteam/vg) and [rpvg](https://github.com/jonassibbesen/rpvg). The spliced pangenome graphs and pantranscriptomes (haplotype-specific transcripts) presented in the paper are avaliable to download in the [Data](#Data) section for use in other projects.

This repository is organized in three subdirectories.

1. The [installation_and_demo](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/installation_and_demo) directory contains installation directions for `vg` and `rpvg`. It also includes a short demo of using the tools for transcriptomic inference, with example data included.

2. The [scripts](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/scripts) directory contains all of the scripts used for analysis and plotting in this project. It is further subdivided by which language the scripts are written in. However, the scripts in the [bash](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/scripts/bash) subdirectory are not the exact scripts we used. They have been simplified to make them easier for others to use, mainly by removing hard-coded paths and replacing environment-defined variables with variables that can be easily edited.

3. The [originals](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/originals) directory contains the raw, unedited bash scripts, as well as the log files. These files are not particularly user-friendly as they include a lot of hard-coded paths. However, we have included them here for transparency and reproducibility. By looking at the scripts and log files you can see exactly how each method was run in the paper. Most of the log files will include a short header which specifies the Docker image that was used. The Docker containers are available at this repository: [s3script-dockerfiles](https://github.com/jonassibbesen/s3script-dockerfiles). For the log files without this header it should be clear from the script itself what version was used. 

4. The [dockerfiles](https://github.com/jonassibbesen/vgrna-project-paper/tree/main/dockerfiles) directory contains recommended Docker files for running scripts in this repository.

## Data

Here you can find links to the data used in the paper. This includes both raw data and data constructed as part of the analyses in the paper. The constructed data included here is data that are either not guaranteed to be reproducible (subsampled transcript annotation and simulated reads) or that are deemed potentially useful in other projects (graphs, pantranscriptomes and indexes).

### Graphs, pantranscriptomes and indexes

The spliced pangenome graphs, pantranscriptomes and indexes:

* http://cgl.gi.ucsc.edu/data/vgrna/pantranscriptomes/

### Genome

The GRCh38 (primary assembly) reference genome:

* ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/

### Transcripts

The GENCODE v29 (primary assembly) transcript annotation:

* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz 

The subsampled (80%) GENCODE v29 transcript annotation:

* http://cgl.gi.ucsc.edu/data/vgrna/transcript_annotation/

### Variants and haplotypes

The 1000 Genomes Project variants and haplotypes lifted to GRCh38:

* http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ 

The IPD-IMGT/HLA gene allele sequences: 

* https://github.com/ANHIG/IMGTHLA (release 3.43.0)

### Reads

The simulated RNA-seq reads:

* http://cgl.gi.ucsc.edu/data/vgrna/simulated_data/

The real RNA-seq reads:

* https://www.ncbi.nlm.nih.gov/sra/?term=SRR1153470 (NA12878)
* https://www.encodeproject.org/experiments/ENCSR000AED/ (NA12878, replicate 1)
* https://github.com/nanopore-wgs-consortium/CHM13/ (CHM13, replicate 1)
* https://www.encodeproject.org/experiments/ENCSR146ZKR/ (Adrenal gland)
* https://www.encodeproject.org/experiments/ENCSR825GWD/ (Sigmoid colon)
* https://www.encodeproject.org/experiments/ENCSR686JJB/ (Adipose tissue)
* https://www.encodeproject.org/experiments/ENCSR502OTI/ (Psoas muscle)
* https://www.encodeproject.org/experiments/ENCSR995BHD/ (Aorta)
* https://www.internationalgenome.org/data-portal/data-collection/geuvadis (NA07051, NA11832, NA11840, NA11930, NA12287, NA12775, NA06994, NA07037, NA07357, NA11829, NA11893, NA12006, NA12043, NA12234, NA12272 & NA12275)
* https://www.internationalgenome.org/data-portal/data-collection/structural-variation (NA19238, NA19239, NA19240, HG00512, HG00513, HG00514, HG00731, HG00732 & HG00733)
* https://www.ncbi.nlm.nih.gov/sra/?term=SRR12765650 (intestine, African American individual)
* https://www.ncbi.nlm.nih.gov/sra/?term=SRR12765534 (intestine, European American individual)

The Iso-Seq alignments:

* https://www.encodeproject.org/experiments/ENCSR706ANY/ (NA12878, all replicates)

