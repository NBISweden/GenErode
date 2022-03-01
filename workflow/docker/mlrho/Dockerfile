FROM ubuntu:18.04 AS builder
MAINTAINER NBIS Sweden <generode@nbis.se>
RUN apt update -y && apt install -y \
	autoconf \
	automake \
	gcc \
	libbz2-dev \
	libcurl4-gnutls-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libssl-dev \
	make \
	perl \
	wget \
	zlib1g-dev

WORKDIR /mlRho
RUN wget http://guanine.evolbio.mpg.de/mlRho/mlRho_2.9.tgz && \
	tar xfvz mlRho_2.9.tgz && \
	cd MlRho_2.9 && \
	make && \
	cp mlRho /mlRho && \
	cd /mlRho

RUN wget http://guanine.evolbio.mpg.de/mlRho/formatPro_0.5.tgz && \
	tar xfvz formatPro_0.5.tgz && \
	cd FormatPro_0.5 && \
	make && \
	cp formatPro /mlRho && \
	cd /mlRho

RUN wget http://guanine.evolbio.mpg.de/mlRho/sam2pro_0.8.tgz && \
	tar xvfz sam2pro_0.8.tgz && \
	cd Sam2pro_0.8 && \
	make && \
	cp sam2pro /mlRho && \
	cd /mlRho

###################################################################################################################

FROM ubuntu:18.04 AS production
RUN apt update -y && apt install -y \
	libcurl4-gnutls-dev \
	libgsl0-dev
COPY --from=nbisweden/generode-samtools_1.9:latest /samtools/samtools /usr/bin
COPY --from=builder /mlRho/mlRho /usr/bin/
COPY --from=builder /mlRho/formatPro /usr/bin/
COPY --from=builder /mlRho/sam2pro /usr/bin
