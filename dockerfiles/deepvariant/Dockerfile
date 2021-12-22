FROM google/deepvariant:1.2.0

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install essential tools

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	wget time git unzip && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip
RUN	pip3 install --no-cache-dir requests awscli

### Add s3 run script

ADD run.sh /home
