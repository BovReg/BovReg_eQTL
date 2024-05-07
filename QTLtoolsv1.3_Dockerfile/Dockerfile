
###############################################################
# Dockerfile to build container images for QTLtools (v1.3.1)
# Based on ubuntu
###############################################################

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies first
RUN apt-get -qq update && \
    apt-get -qq -y install \
    build-essential gcc g++ make libtool texinfo dpkg-dev pkg-config \
    libgsl-dev wget locate less vim zlib1g-dev bzip2 lzma curl r-base \
    libboost-dev libcurl4-openssl-dev libboost-all-dev libbz2-dev liblzma-dev \
    libpcre3 libpcre3-dev

# Install vcftools and bcftools

RUN  apt-get  update
RUN  apt-get -y  install gawk 

RUN  apt-get -y  install vcftools

RUN  apt-get -y  install bcftools 


# Install R (v.3.6.3)
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz && \
    tar -zxvf R-3.6.3.tar.gz && \
    rm -rf R-3.6.3.tar.gz && \
    cd R-3.6.3 && \
    ./configure --with-x=no --with-readline=no && \
    cd src/nmath/standalone/ && \
    make

# Install R library
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install('qvalue')"

# Install HTSlib (v.1.3.1)
RUN wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar xjvf htslib-1.3.1.tar.bz2 && \
    rm -rf htslib-1.3.1.tar.bz2 && \
    cd /htslib-1.3.1/ && \
    ./configure && \
    make

# Install QTLtools (v1.3.1)
RUN wget https://qtltools.github.io/qtltools/binaries/QTLtools_1.3.1_source.tar.gz && \
    tar -zxvf QTLtools_1.3.1_source.tar.gz && \
    rm -rf QTLtools_1.3.1_source.tar.gz && \
    cd /qtltools/ && \
    sed \
    -e 's/BOOST_INC=/BOOST_INC=\/usr\/include/' \
    -e 's/BOOST_LIB=/BOOST_LIB=\/usr\/lib\/x86_64-linux-gnu/' \
    -e 's/RMATH_INC=/RMATH_INC=\/R-3.6.3\/src\/include/' \
    -e 's/RMATH_LIB=/RMATH_LIB=\/R-3.6.3\/src\/nmath\/standalone/' \
    -e 's/HTSLD_INC=/HTSLD_INC=\/htslib-1.3.1/' \
    -e 's/HTSLD_LIB=/HTSLD_LIB=\/htslib-1.3.1/' \
    /qtltools/Makefile > /qtltools/Makefile_2 && \
    rm -rf /qtltools/Makefile && \
    mv /qtltools/Makefile_2 /qtltools/Makefile && \
    make && \
    make install

# Install bedtools

RUN apt-get install bedtools






# Set working directory
WORKDIR /home

