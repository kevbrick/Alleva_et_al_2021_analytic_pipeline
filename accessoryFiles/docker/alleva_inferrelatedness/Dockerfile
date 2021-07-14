FROM nfcore/base:1.10.2

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for inferring PRDM9 template switches"

## Install gcc
RUN apt-get update
RUN apt-get install -y make perl cpanminus gcc

RUN cpanm -v List::Util
RUN cpanm -v List::MoreUtils
RUN cpanm -v List::UtilsBy
RUN cpanm -v Getopt::Long
RUN cpanm -v String::LCSS
RUN cpanm -v String::LCSS_XS

ENV PERL5LIB /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2:$PERL5LIB

RUN wget https://cpan.metacpan.org/authors/id/G/GR/GROMMEL/Math-Round-0.07.tar.gz -P /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2
RUN tar -zxvf /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math-Round-0.07.tar.gz
RUN mv Math-Round-0.07 /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

