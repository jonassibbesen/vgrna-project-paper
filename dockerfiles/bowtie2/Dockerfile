FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install bowtie2 

RUN wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-linux-x86_64.zip/download && \
	unzip download && \
	cd bowtie2-2.4.1-linux-x86_64 && \
	mv bowtie2 /usr/bin/ && \
	mv bowtie2-* /usr/bin/ && \
	mkdir -p /scripts/bowtie2/ && \
	mv scripts/* /scripts/bowtie2/ && \
	cd .. && \
	rm -r download* && \
	rm -r bowtie2-2.4.1*
