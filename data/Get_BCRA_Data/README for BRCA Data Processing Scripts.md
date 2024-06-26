# README for BRCA Data Processing Scripts

## Overview

This repository contains scripts for processing and obtaining methylation data relevant to breast cancer research. The primary outputs of these scripts are two files: `BRCA_data_meth.csv` ,`BRCA_data_meth_range.csv`, `need_gene_exp_13988(13982).csv`. These files are generated by following methodologies adapted from the study "Collective effects of long-range DNA methylations predict gene expressions and estimate phenotypes in cancer".

## Files Description

- **BRCA_data_meth.csv**: This file contains methylation data for breast cancer. The methylation data represents epigenetic modifications that can influence gene expression without altering the underlying DNA sequence.
- **BRCA_data_meth_range.csv**: This file provides information on the range of methylation sites. It details the genomic coordinates of methylation across the genome, which can be critical for understanding the locational impact of methylation on gene regulation.
- **need_gene_exp_13988(13982).csv**: This file contains gene expression data for breast cancer. Six of the genes are duplicated.

## Data Acquisition

To generate these files, the script `function.study.name.download_files` is executed, which:

1. Downloads necessary datasets from specified URLs.
2. Processes the data to extract relevant methylation and genomic range information.
3. Saves the processed data into CSV files for further analysis.

## Usage

To run the scripts and generate the required `.csv` files, follow these steps:

1. Ensure that all required libraries, such as `GenomicRanges` and `glmnet`, are installed in your R environment.
2. Run the script `function.study.name.download_files` by sourcing it in your R console or script environment.
3. The scripts will automatically process the downloaded data and save them in the current working directory as  `BRCA_data_meth.csv` ,`BRCA_data_meth_range.csv`, `need_gene_exp_13988(13982).csv`.

## Dependencies

- R (version 3.6 or higher recommended)
- GenomicRanges: For handling genomic intervals and coordinates.
- glmnet: For modeling data using regularized regression.

Please ensure that your R environment is set up with the required libraries before executing the scripts. 

