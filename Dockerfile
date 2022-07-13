FROM rocker/tidyverse:4.2

LABEL maintainer="Daniel Miller (millerd15@chop.edu)"

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev curl

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk

# Install required R packages from CRAN
RUN install2.r cutpointr \
	ggforce \
	openxlsx \
	patchwork \
	R.utils \
	survminer \
	Hmisc

# Install maftools from BioC
RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r maftools

# Install R packages from GitHub
RUN installGithub.r d3b-center/annoFuse \
	jokergoo/ComplexHeatmap

ADD Dockerfile .
