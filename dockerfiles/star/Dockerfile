FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install star 

RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
	tar -xvzf 2.7.9a.tar.gz && \
	cd STAR-2.7.9a && \
	mv bin/Linux_x86_64_static/STAR* /usr/bin/ && \
	mkdir -p /scripts/star/ && \
	mv extras/scripts/* /scripts/star/ && \
	cd .. && \
 	rm -r *2.7.9a*
