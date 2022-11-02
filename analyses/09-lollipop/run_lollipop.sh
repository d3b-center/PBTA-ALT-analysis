#!/bin/bash

set -e
set -o pipefail

# OncoKB access token required to run module; if not provided then exit
if [ "$ONCOKB" = "" ]; then
    echo 'Error: An Oncokb access token is required to run this module in full. Please provide token before running analysis module in the format ONCOKB=<TOKEN> bash run_lollipop.sh'
    exit 1
fi

# create subset MAF
Rscript -e "rmarkdown::render('code/01-create-subset-maf.Rmd')"

# Run oncokb MAF annotation: 
ONCOKB=$ONCOKB bash code/02-run-oncokb-annotation.sh

# make lollipop plots
Rscript -e "rmarkdown::render('code/03-create-lollipop.Rmd')"
