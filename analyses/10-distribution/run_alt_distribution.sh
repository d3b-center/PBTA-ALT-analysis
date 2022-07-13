#!/bin/bash

set -e
set -o pipefail

# ALT distribution plots
Rscript 01-alt-distribution.R
