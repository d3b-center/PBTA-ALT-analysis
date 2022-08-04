#!/bin/bash

set -e
set -o pipefail

# ALT distribution plots
Rscript code/00-table2-somatic-mmr.R
Rscript code/01-tables.R

