FROM ubuntu:18.04 AS builder
MAINTAINER NBIS Sweden <generode@nbis.se>
RUN apt update -y && apt install -y \
	autoconf \
	automake \
	gcc \
	libbz2-dev \
	libcurl4-gnutls-dev \
	liblzma-dev \
	libncurses5-dev \
	libssl-dev \
	make \ 
	perl \
	wget \ 
	zlib1g-dev

WORKDIR /htslib
RUN 	wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 && \
	bunzip2 htslib-1.15.1.tar.bz2 && \
	tar xfv htslib-1.15.1.tar && \
	cd htslib-1.15.1 && \
	./configure && \
	make && \
	cp bgzip /htslib/ && \
	cp tabix /htslib/ && \
	make install
