#!/bin/bash

set -e
set -o pipefail

# get dna alterations
Rscript -e "rmarkdown::render('code/01_get_alt_TERT_ATRX_DAXX.Rmd')"

# get rnaseq expression and fusion
Rscript code/02_get_alt_rnaseq.R

# merge with telomere content ,ccircle and terra and plot
Rscript -e "rmarkdown::render('code/03_merge_telomere_content.Rmd')"
