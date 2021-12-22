FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install salmon 

RUN wget --no-check-certificate https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz && \
	tar -xvzf salmon-1.5.2_linux_x86_64.tar.gz && \
	cd salmon-1.5.2_linux_x86_64 && \
	mv bin/salmon /usr/bin/ && \
	mv lib/* /usr/lib/ && \
	cd .. && \
 	rm -r salmon*
