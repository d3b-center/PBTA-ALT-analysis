#!/bin/bash

set -e
set -o pipefail

# ALT distribution plots
Rscript code/01-alt-distribution.R
