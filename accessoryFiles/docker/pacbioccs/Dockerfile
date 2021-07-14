FROM nfcore/base:1.10.2

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for calling PRDM9 haplotypes from long read data"

## Install gcc
RUN apt-get update
RUN apt-get install -y make

COPY environment.yml /
RUN conda env create --quiet -f environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/pbtools/bin:$PATH

