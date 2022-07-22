#!/bin/bash

set -e
set -o pipefail

# create oncoprint matrix
Rscript code/00-create-subset-matrix.R

# create oncoprint
Rscript code/01-plot_oncoplot.R
