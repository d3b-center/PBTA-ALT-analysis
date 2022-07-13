#!/bin/bash

set -e
set -o pipefail

# run cutpoint analysis
Rscript -e "rmarkdown::render('code/cutpoint-analysis.Rmd')"