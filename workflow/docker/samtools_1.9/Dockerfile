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
WORKDIR /samtools
RUN 	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	bunzip2 samtools-1.9.tar.bz2 && \
	tar xfv samtools-1.9.tar && \
	cd samtools-1.9 && \
	./configure && \
	make && \
	cp samtools /samtools/ && \
	make install
