# Set the base image to Ubuntu 16.04 and NVIDIA GPU
FROM nvidia/cuda:11.1-cudnn8-devel-ubuntu16.04

# File Author / Maintainer
MAINTAINER Kevin Brick <kevbrick@gmail.com>

ARG PACKAGE_VERSION=4.0.15
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes $BUILD_PACKAGES \
                          libcurl4-openssl-dev \
                          libssl-dev \
                          libhdf5-cpp-11 \
                          libzmq5 \
                          libboost-atomic1.58.0 \
                          libboost-chrono1.58.0 \
                          libboost-date-time1.58.0 \
                          libboost-filesystem1.58.0 \
                          libboost-program-options1.58.0 \
                          libboost-regex1.58.0 \
                          libboost-system1.58.0 \
                          libboost-log1.58.0 \
                          libboost-iostreams1.58.0 \
                          wget && \
    cd /tmp &&\
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${PACKAGE_VERSION}-1~xenial_amd64.deb && \
    dpkg -i --ignore-depends=nvidia-384,libcuda1-384 /tmp/ont_guppy_${PACKAGE_VERSION}-1~xenial_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    curl \
    git \
    python3.5 \
    python3-pip

RUN python3 -m pip install --upgrade pip

RUN pip install ont-fast5-api==3.1.6

RUN pip install ont-bonito==0.2.3




