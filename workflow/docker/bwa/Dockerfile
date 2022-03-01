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

WORKDIR /bwa
RUN	wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
        bunzip2 bwa-0.7.17.tar.bz2 && \
	tar xfv bwa-0.7.17.tar && \
	cd bwa-0.7.17 && \
	make && \
	cp bwa /bwa/

FROM ubuntu:18.04 AS production
RUN apt update -y && apt install -y \
	libcurl4-gnutls-dev
COPY --from=nbisweden/generode-samtools_1.9:latest /samtools/samtools /usr/bin
COPY --from=builder /bwa/bwa /usr/bin/
