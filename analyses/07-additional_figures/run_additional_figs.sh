#!/bin/bash

set -e
set -o pipefail

# investigate subtype distributions
Rscript code/mol_subtype_bar.R

# mutation counts by gene
Rscript code/mut_count_box.R

# mutation counts frequencies
Rscript code/mut_count_freq.R