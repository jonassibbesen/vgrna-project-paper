FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install rsem 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar -xvzf v1.3.3.tar.gz && \
    cd RSEM-1.3.3 && \
    make install && \
    mkdir -p /scripts/rsem/ && \
    mv EBSeq /scripts/rsem/ && \
    mv pRSEM /scripts/rsem/ && \
    cd .. && \
    rm -r *1.3.3*

### Install bowtie2 

RUN wget --no-check-certificate https://github.com/BenLangmead/bowtie2/releases/download/v2.4.4/bowtie2-2.4.4-linux-x86_64.zip && \
    unzip bowtie2-2.4.4-linux-x86_64.zip && \
    cd bowtie2-2.4.4-linux-x86_64 && \
    mv bowtie2 /usr/bin/ && \
    mv bowtie2-* /usr/bin/ && \
    mkdir -p /scripts/bowtie2/ && \
    mv scripts/* /scripts/bowtie2/ && \
    cd .. && \
    rm -r bowtie2-2.4.4-linux-x86_64*
