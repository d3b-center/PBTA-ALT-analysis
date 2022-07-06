FROM rocker/tidyverse:4.2

LABEL maintainer="Daniel Miller (millerd15@chop.edu)"

RUN apt update && apt install -y zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev

RUN install2.r cutpointr ggforce openxlsx patchwork R.utils survminer

RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r maftools

RUN installGithub.r d3b-center/annoFuse jokergoo/ComplexHeatmap

ADD Dockerfile .
