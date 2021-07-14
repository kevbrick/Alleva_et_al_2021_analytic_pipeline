FROM r-base

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for gel quantification script"

## Install gcc
RUN apt-get update
RUN apt-get install -y make
RUN apt-get install -y perl
RUN apt-get -y install libblas-dev bioperl 
RUN apt-get install -y cpanminus procps
RUN apt-get install -y libfftw3-dev libfftw3-3 libtiff-dev pngtools 
RUN apt-get install -y imagemagick libmagic-dev libmagick++-6.q16-dev libmagickcore-6.q16-dev libmagickcore-dev libmagick++-dev libmagickwand-6.q16-dev libmagickwand-dev libmagic-ocaml-dev libmagics++-dev libmagics++-metview-dev

RUN cpanm -v List::Util
RUN cpanm -v List::MoreUtils
RUN cpanm -v List::UtilsBy
RUN cpanm -v Getopt::Long
RUN cpanm -v String::LCSS

ENV PERL5LIB /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2:$PERL5LIB

RUN wget https://cpan.metacpan.org/authors/id/G/GR/GROMMEL/Math-Round-0.07.tar.gz -P /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2
RUN tar -zxvf /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math-Round-0.07.tar.gz
RUN mv Math-Round-0.07 /opt/conda/envs/parsePRDM9RE/lib/site_perl/5.26.2/Math

RUN R -e "if('imager'        %in% rownames(installed.packages()) == FALSE){install.packages('imager',dependencies=TRUE, repos='http://cran.rstudio.com/')}"

RUN R -e "if('data.table'   %in% rownames(installed.packages()) == FALSE){install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('dplyr'        %in% rownames(installed.packages()) == FALSE){install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('extrafont'    %in% rownames(installed.packages()) == FALSE){install.packages('extrafont',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('ggplot2'      %in% rownames(installed.packages()) == FALSE){install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('ggpmisc'      %in% rownames(installed.packages()) == FALSE){install.packages('ggpmisc',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('ggpubr'       %in% rownames(installed.packages()) == FALSE){install.packages('ggpubr',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('ggrepel'      %in% rownames(installed.packages()) == FALSE){install.packages('ggrepel',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('grid'         %in% rownames(installed.packages()) == FALSE){install.packages('grid',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('gridExtra'    %in% rownames(installed.packages()) == FALSE){install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')}"

RUN R -e "if('lsr'          %in% rownames(installed.packages()) == FALSE){install.packages('lsr',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('numform'      %in% rownames(installed.packages()) == FALSE){install.packages('numform',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('plyr'         %in% rownames(installed.packages()) == FALSE){install.packages('plyr',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('png'          %in% rownames(installed.packages()) == FALSE){install.packages('png',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('reshape2'     %in% rownames(installed.packages()) == FALSE){install.packages('reshape2',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('stringi'      %in% rownames(installed.packages()) == FALSE){install.packages('stringi',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('stringr'      %in% rownames(installed.packages()) == FALSE){install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')}"
RUN R -e "if('tidyverse'    %in% rownames(installed.packages()) == FALSE){install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')}"

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

