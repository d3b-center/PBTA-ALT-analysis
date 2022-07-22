#!/bin/bash

set -e
set -o pipefail

# get pedcbio ids for comparisons
Rscript code/gather-pedc-ids.R