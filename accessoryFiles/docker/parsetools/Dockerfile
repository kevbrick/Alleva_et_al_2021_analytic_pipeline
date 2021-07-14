FROM python:3.7

LABEL authors="Kevin Brick" \
      description="Docker image containing parsing tools"

## Install gcc
RUN apt-get update && \
    apt-get -y --force-yes install make wget && \
    python3 -m pip install xlsx2csv && \
    python3 -m pip install pdf2text
