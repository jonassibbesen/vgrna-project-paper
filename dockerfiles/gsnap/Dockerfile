FROM jsibbesen/base-s3script:18.04-s1

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install gsnap 

RUN wget --no-check-certificate http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-06-30.tar.gz && \
	tar -xvzf gmap-gsnap-2020-06-30.tar.gz && \
	cd gmap-2020-06-30 && \
	./configure && \
	make && \
	make install && \
	cd .. && \
 	rm -rf gmap-*
