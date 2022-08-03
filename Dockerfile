FROM rocker/tidyverse:4.2

LABEL maintainer="Jo Lynne Rokita (rokita@chop.edu)"

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev curl

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk \
   libxt6

# Install required R packages from CRAN
RUN install2.r \
	ggforce \
	openxlsx \
	patchwork \
	R.utils \
	survminer \
	Hmisc \
	optparse \
	cutpointr

# Install R packages from GitHub
RUN installGithub.r d3b-center/annoFuse \
	jokergoo/ComplexHeatmap \ 
	PoisonAlien/maftools

# Install pip3 and python reqs for oncokb
RUN apt-get -y --no-install-recommends install \
    python3-pip python3-dev
RUN pip3 install \
  "matplotlib==3.1.2" \
  "kiwisolver==1.2.0"

# Install oncokb
RUN git clone https://github.com/oncokb/oncokb-annotator.git /home/oncokb-annotator

ADD Dockerfile .
