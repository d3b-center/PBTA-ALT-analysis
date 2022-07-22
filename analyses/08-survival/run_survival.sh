#!/bin/bash

set -e
set -o pipefail

# survival analysis
Rscript -e "rmarkdown::render('code/01-survival-analysis_subtypes.Rmd')"

# kaplan meier
Rscript code/02-kaplan-meier.R

# forest plot
Rscript code/03-hgg-alt-subtype-forest-plot.R