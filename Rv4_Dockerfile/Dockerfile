FROM rocker/r-ubuntu:18.04

MAINTAINER Joachim Gassen "gassen@wiwi.hu-berlin.de" 
# Update packages lists

RUN apt-get update 
# Latex for pdflatex (needed for R CMD check)

RUN apt-get install -y --no-install-recommends\
  texlive-latex-base texlive-fonts-recommended texlive-fonts-extra \
  texlive-latex-extra texinfo
# Extra stuff that your package or its dependencies need

RUN apt-get install -y --no-install-recommends \
  libxml2-dev libcurl4-openssl-dev libssl-dev
# Update R to 4.0 - from
# http://dirk.eddelbuettel.com/papers/r4_r400_binaries.pdf
# Update to current repo state

RUN apt-get dist-upgrade -y
# then add the edd/r-4.0 PPA

RUN add-apt-repository -y ppa:edd/r-4.0
# upgrade again

RUN apt-get dist-upgrade -y
# Add the R remotes package to eventually install your package if it resides
# on Github

RUN  R -e "install.packages('remotes')"
# Install dependencies of your package (yes.. {tidycovid19} _is_
# dependency heavy ...)

RUN R -e "install.packages(c('xml2', 'rvest', 'readxl', \
 'magrittr', 'dplyr', 'tidyr', 'countrycode', \
 'ggplot2', 'wbstats', 'gghighlight', 'gtrendsR', \
 'ggrepel', 'readr', 'rlang', 'scales', 'lubridate', \
 'shiny', 'shinyWidgets', 'shinyjs', 'rclipboard', \
 'stringr', 'zoo', 'maps', 'gganimate'))"

 


FROM bioconductor/bioconductor_docker:devel



LABEL authors="Praveen" \
      description="Docker image containing all requirements for the eQTLUtils R package and pipeline"

RUN R -e "BiocManager::install(c('BiocCheck','SummarizedExperiment','lumi', 'limma', 'dplyr','cqn','ggplot2','htmlwidgets', 'tidyr','assertthat','devtools','GenomicRanges','readr', 'GDSArray','plotly','SNPRelate'))"








