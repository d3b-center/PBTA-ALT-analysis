# ALT in Pediatric Brain Tumors Can Occur without ATRX Mutation and is Enriched in Patients with Pathogenic Germline MMR Variants
Jennifer L. Stundon, Heba Ijaz, Krutika S. Gaonkar, Rebecca S. Kaufman, Run Jin, Anastasios Karras, Zalman Vaksman, Jung Kim, Daniel P. Miller, Yiran Guo, Mariarita Santi, Marilyn Li, Gonzalo Lopez, Phillip J. Storm Adam C. Resnick, Angela Waanders, Suzanne McFarland, Douglas R. Stewart, Sharon J. Diskin, Jo Lynne Rokita\*, Kristina A. Cole
\*Corresponding author


## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct somatic analyses for the manuscript noted above.

1. Clone the repository
```
git clone https://github.com/d3b-center/PBTA-ALT-analysis.git
```

2. Pull the docker container:
```
docker pull jrokita1/pbta-alt:version1.0
```

3. Start the docker container, from the `PBTA-ALT-analysis` folder, run:
```
docker run --name container_name -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/PBTA-ALT-analysis pbta-alt:latest
```

4. To execute shell within the docker image, from the `PBTA-ALT-analysis` folder, run:
```
docker exec -ti container_name bash
```

5. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/PBTA-ALT-analysis/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
│   ├── 01-alterations_ratio_check
│   ├── 02-add-histologies
│   ├── 03-get-pedcbio-ids
│   ├── 04-cutpoint-analysis
│   ├── 05-oncoplot
│   ├── 07-additional_figures
│   ├── 08-survival
│   ├── 09-lollipop
│   ├── 10-distribution
│   └── 11-tables
├── data
│   ├── TableS2-DNA-results-table.xlsx -> v1/TableS2-DNA-results-table.xlsx
│   ├── TableS3-RNA-results-table.xlsx -> v1/TableS3-RNA-results-table.xlsx
│   ├── consensus_wgs_plus_cnvkit_wxs.tsv.gz -> v1/consensus_wgs_plus_cnvkit_wxs.tsv.gz
│   ├── histologies.tsv -> v1/histologies.tsv
│   ├── pbta-fusion-putative-oncogenic.tsv -> v1/pbta-fusion-putative-oncogenic.tsv
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds -> v1/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds -> v1/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
│   ├── pbta-histologies.tsv -> v1/pbta-histologies.tsv
│   ├── pbta-snv-consensus-mutation-tmb-coding.tsv -> v1/pbta-snv-consensus-mutation-tmb-coding.tsv
│   ├── pbta-sv-manta.tsv.gz -> v1/pbta-sv-manta.tsv.gz
│   ├── release-notes.md -> v1/release-notes.md
│   ├── snv-consensus-plus-hotspots.maf.tsv.gz -> v1/snv-consensus-plus-hotspots.maf.tsv.gz
│   └── v1
├── download-data.sh
└── scratch
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)), Krutika S. Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Run Jin ([@runjin326](https://github.com/runjin326)), and Daniel P. Miller ([@dmiller15](https://github.com/dmiller15))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita: rokita@chop.edu

