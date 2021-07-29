FROM python:3.8

LABEL authors="Kevin Brick" \
      description="Docker image containing parsing tools"

## Install gcc
RUN apt-get update 
##   apt-get -y --force-yes install make wget

RUN python3 -m pip install pdf2text

RUN python3 -m pip install xlsx2csv
    

RUN apt-get -y install tesseract-ocr
