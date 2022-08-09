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

WORKDIR /bedtools

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
        mv bedtools.static.binary /bedtools/bedtools && \
        chmod a+x /bedtools/bedtools

FROM ubuntu:18.04 AS production
RUN apt update -y && apt install -y \
	libcurl4-gnutls-dev
# COPY --from=nbisweden/generode-htslib-1.15.1:latest /htslib/bgzip /usr/bin
# COPY --from=nbisweden/generode-htslib-1.15.1:latest /htslib/tabix /usr/bin
COPY --from=verku/htslib-1.15.1:latest /htslib/bgzip /usr/bin
COPY --from=verku/htslib-1.15.1:latest /htslib/tabix /usr/bin
COPY --from=builder /bedtools/bedtools /usr/bin/
