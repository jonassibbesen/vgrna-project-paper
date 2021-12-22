FROM quay.io/vgteam/vg:ci-42-385fd63622288dbeb4da76f83b16667bb350f397

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install essential tools

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	wget time git unzip autoconf automake make cmake gcc g++ gdb \
	python2.7 python-setuptools python3 python3-pip python3-setuptools && \
    rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python2 /usr/bin/python

RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir requests awscli

### Install htslib 

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev && \
    rm -rf /var/lib/apt/lists/*
    
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
	tar -xvjf htslib-1.9.tar.bz2 && \
	cd htslib-1.9 && \
	./configure && \
	make install && \
 	cd .. && \
 	rm -r *

### Install bcftools 

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	libperl-dev libgsl0-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
	tar -xvjf bcftools-1.9.tar.bz2 && \
	cd bcftools-1.9 && \
	./configure && \
	make install && \
 	cd .. && \
 	rm -r *

### Install samtools 

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	libncurses5-dev && \
    rm -rf /var/lib/apt/lists/*
    
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar -xvjf samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure && \
	make install && \
 	cd .. && \
 	rm -r *

### Install seqtk 
    
RUN wget --no-check-certificate https://github.com/lh3/seqtk/archive/v1.3.tar.gz && \
	tar -xvzf v1.3.tar.gz && \
	cd seqtk-1.3 && \
	make install && \
 	cd .. && \
 	rm -r *

### Install bedtools 
    
RUN wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
	tar -xvzf bedtools-2.29.1.tar.gz && \
	cd bedtools2 && \
	make && \
	mv bin/* /usr/bin/ && \
	mkdir -p /scripts/bedtools/ && \
	mv scripts/* /scripts/bedtools/ && \
 	cd .. && \
 	rm -r *

### Add s3 run script

ADD run.sh /home
