FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install hisat2 

RUN wget --no-check-certificate https://cloud.biohpc.swmed.edu/index.php/s/4pMgDq4oAF9QCfA/download && \
	unzip download && \
	cd hisat2-2.2.1 && \
	mv hisat2 /usr/bin/ && \
	mv hisat2-* /usr/bin/ && \
	mv *py /usr/bin/ && \
	mkdir -p /scripts/hisat2/ && \
	mv scripts/* /scripts/hisat2/ && \
	cd .. && \
	rm -r download* && \
	rm -r hisat2-2.2.1*
