#!/bin/bash

set -e
set -o pipefail

# lollipop plots
Rscript -e "rmarkdown::render('code/01-create-lollipop.Rmd')"
