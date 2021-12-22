FROM quay.io/jsibbesen/vgdev-s3script:vgdev-385fd636-s2

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install R

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

### Install pysam and biopython

RUN pip3 install --upgrade pip
RUN	pip3 install --no-cache-dir pysam biopython

### Install vgrna scripts 

RUN git clone https://github.com/jonassibbesen/vgrna-project-scripts.git && \
	cd vgrna-project-scripts && \
	git checkout 4b55c11079907761923e8dc41e06f01fb2b91fc6 && \
	cd c++ && \
	cmake . && \
	make && \
	mv bin/* /usr/bin/ && \
	cd .. && \
	mkdir -p /scripts/vgrna/ && \
	mv python/*.py /scripts/vgrna/ && \
	mv R/*.R /scripts/vgrna/ && \
	mv .git/ /scripts/vgrna/ && \
	cd .. && \
 	rm -r vgrna-project-scripts
