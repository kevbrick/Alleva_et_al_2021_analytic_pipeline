FROM nfcore/base:1.10.2

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for associating SNPs with PRDM9 alleles"

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

RUN cd /usr/local && \
  wget https://github.com/vcftools/vcftools/archive/refs/tags/v0.1.16.tar.gz && \
  tar -zxvf v0.1.16.tar.gz && \
  cd /usr/local/vcftools-0.1.16 && \
  bash ./autogen.sh && \
  bash ./configure --prefix=/usr/local/vcftools-0.1.16 && \
  make && \
  make install && \
  make clean
  
RUN cd /usr/local && \
  wget http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip && \
  unzip plink-1.07-x86_64.zip 

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

RUN apt-get install -y gcc python3-dev python3-setuptools 
RUN python3 -m pip install --no-cache-dir -U crcmod

RUN cd /usr/local && \
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-347.0.0-linux-x86_64.tar.gz && \
    tar -zxvf google-cloud-sdk-347.0.0-linux-x86_64.tar.gz && \
    ./google-cloud-sdk/install.sh

ENV PATH=/usr/local/bedtools2/bin:/usr/local/vcftools-0.1.16/bin:/usr/local/plink-1.07-x86_64:/usr/local/samtools/bin:/usr/local/htslib/bin:$PATH
