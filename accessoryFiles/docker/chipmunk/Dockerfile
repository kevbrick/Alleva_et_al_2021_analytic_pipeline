FROM ruby:2.5.8-alpine3.12

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for chipmunk motif calling"

RUN apk add openjdk8-jre
RUN apk add imagemagick
RUN apk add imagemagick-dev
RUN apk add g++
RUN apk add make
RUN apk add bash
RUN apk add perl

## Install gcc
#RUN apt-get update
#RUN apt-get -y --force-yes install libmagickwand-dev imagemagick ruby-full

RUN gem install rmagick

RUN mkdir /opt/chipmunk
RUN mkdir /opt/sarus
RUN wget -O /opt/chipmunk/chipmunk_scripts_v8upd.zip http://autosome.ru/ChIPMunk/chipmunk_scripts_v8upd.zip
RUN wget -O /opt/chipmunk/userguide_2020.pdf         http://autosome.ru/ChIPMunk/userguide_2020.pdf
RUN wget -O /opt/sarus/sarus.jar                     http://autosome.ru/ChIPMunk/sarus.jar

RUN unzip /opt/chipmunk/chipmunk_scripts_v8upd.zip -d /opt/chipmunk/

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/chipmunk:/opt/sarus:$PATH
ENV CHIPMUNK=/opt/chipmunk
ENV SARUSJAR=/opt/sarus/sarus.jar

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
