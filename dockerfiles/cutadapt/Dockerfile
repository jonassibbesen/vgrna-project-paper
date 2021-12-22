FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install cutadapt

RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir cutadapt

### Install miRNA adapter scripts

RUN git clone https://github.com/sorenar/miRNA-seq-adapters.git && \
	cd miRNA-seq-adapters && \
	git checkout 6cc7e3833a5e5416908f9b2dc0f76bbe845e0e1d && \
	mkdir -p /scripts/mirna/ && \
	mv * /scripts/mirna/ && \
	cd .. && \
 	rm -r miRNA-seq-adapters
