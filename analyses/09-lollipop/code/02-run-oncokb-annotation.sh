#!/bin/bash

set -e
set -o pipefail

python_path=(/usr/bin/python3)
maf_annotator=(/home/oncokb-annotator/MafAnnotator.py)
dgd_in=(output/snv-dgd-goi.maf.tsv)
hgat_in=(output/hgat-goi.maf.tsv)
dgd_oncokb_out=(output/dgd_maf-goi-oncokb.tsv)
hgat_oncokb_out=(output/snv-consensus-plus-hotspots-hgat-oncokb.maf.tsv)

#OncoKB access token required; if not provided then exit 
if [ "$ONCOKB" = "" ]; then
    echo 'Error: An Oncokb access token is required to run this script. Please provide token before running script in the format: ONCOKB=<TOKEN> bash 02-run-oncokb-annotation.sh'
fi

# Run maf_annotator on dgd samples
$python_path $maf_annotator -i $dgd_in -o $dgd_oncokb_out -b $ONCOKB -q hgvsp_short -r GRCh38

#Run maf_annotator on hgat samples
$python_path $maf_annotator -i $hgat_in -o $hgat_oncokb_out -b $ONCOKB -q hgvsp_short -r GRCh38