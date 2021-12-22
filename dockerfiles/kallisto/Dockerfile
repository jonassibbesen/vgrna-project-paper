FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install kallisto 

RUN wget --no-check-certificate https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz && \
	tar -xvzf kallisto_linux-v0.46.2.tar.gz && \
	cd kallisto && \
	mv kallisto /usr/bin/ && \
	cd .. && \
 	rm -r kallisto*
