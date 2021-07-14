FROM nfcore/base:1.10.2

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for validating new PRDM9 alleles"

## Install gcc
RUN apt-get update && \
    apt-get install -y apt-utils && \
    apt-get install -y make perl cpanminus g++ tar zip libghc-zlib-dev autoconf pkg-config && \
    apt-get install -y wget curl libcurl4 libcurl4-openssl-dev

RUN cpanm -v List::Util && \
  cpanm -v List::MoreUtils && \
  cpanm -v List::UtilsBy && \
  cpanm -v Getopt::Long

ENV PERL5LIB /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2:$PERL5LIB

RUN wget https://cpan.metacpan.org/authors/id/G/GR/GROMMEL/Math-Round-0.07.tar.gz -P /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2 && \
  tar -zxvf /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math-Round-0.07.tar.gz && \
  mv Math-Round-0.07 /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math

## GET BEDTOOLS
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
  mkdir -p /usr/local/bedtools2/bin && \
  mv bedtools.static.binary /usr/local/bedtools2/bin/bedtools && \
  echo '#!/bin/bash\n${0%/*}/bedtools slop "$@"' >/usr/local/bedtools2/bin/slopBed && \
  echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/intersectBed && \
  echo '${0%/*}/bedtools intersect "$@"' >>/usr/local/bedtools2/bin/intersectBed && \
  echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/mergeBed && \
  echo '${0%/*}/bedtools merge "$@"'     >>/usr/local/bedtools2/bin/mergeBed && \
  echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/mapBed && \
  echo '${0%/*}/bedtools map "$@"'       >>/usr/local/bedtools2/bin/mapBed && \
  chmod a+x /usr/local/bedtools2/bin/*

RUN apt-get install bzip2 libbz2-dev liblzma-dev && \
  cd /usr/local && \
  wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
  tar -xjvf htslib-1.12.tar.bz2 && \
  cd htslib-1.12 && \
  ./configure --prefix=/usr/local/htslib && \
  make && \
  make install

RUN cd /usr/local && \
  wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
  tar -xjvf samtools-1.12.tar.bz2 && \
  cd samtools-1.12 && \
  ./configure --prefix=/usr/local/samtools && \
  make && \
  make install

RUN cd /usr/local && \
  curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 | tar -jxvf -

ENV PATH=/usr/local/bedtools2/bin:/usr/local/minimap2-2.20_x64-linux:/usr/local/samtools/bin:/usr/local/htslib/bin:$PATH
