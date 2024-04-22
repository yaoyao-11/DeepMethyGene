# README for DeepMethyGene

## Overview

we propose DeepMethyGene, employing conventional deep learning methods like CNN with variable convolutional kernels and ResNet blocks for prediction. Referencing geneEXPLORE, DeepMethyGene was used to predict the expression levels of 13,982 genes in TCGA breast cancer data, achieving a five-fold cross-validation result of 0.64.

## Frame work of DeepMethyGene 

This diagram illustrates the framework of DeepMethyGene. Input Data: The input is gene expression data from 13,892 genes. The number of methylations between the left and right 1Mb of different gene expression sites is different.
<img width="960" alt="DeepMethyGene" src="https://github.com/yaoyao-11/DeepMethyGene/assets/84023156/d07840dc-6f82-47db-a29e-b83f8968e0b3">


## Data Files Description

### need_gene_exp_13988(13982).csv

This file contains gene expression values needed for the analysis. It includes data for 13,988 genes, with 6 genes duplicated within the dataset. To ensure accuracy, the duplicates are resolved by selecting the gene with the most reliable expression results based on the context of the analysis. This selection process is crucial for maintaining the integrity and precision of the gene expression analysis.

### BRCA_data_meth.csv

As previously described, this CSV file includes methylation data specifically related to breast cancer. It compiles information on how methylation patterns are distributed across the genome of breast cancer cells, offering insights into the epigenetic regulation involved in the disease.

### BRCA_data_meth_range.csv

This file details the genomic ranges of methylation sites associated with breast cancer. It is crucial for researchers aiming to understand the spatial distribution of methylation and its potential impact on gene regulation in cancerous cells.

### hg19_promoter.txt

This text file provides data on promoters based on the hg19 human genome assembly. Promoters are critical regions of DNA where transcription of a gene by RNA polymerase begins, hence understanding their locations can help in linking epigenetic modifications like methylation to changes in gene expression.

### gene_list.csv

This CSV file contains a list of specific genes required for the study. The list is used to filter and analyze gene-specific data across different datasets, such as expression levels or methylation patterns, facilitating focused and precise biological insights.

## Usage

To run the scripts and generate the required `.csv` files, follow these steps:

1. Ensure that all required libraries, such as `GenomicRanges` and `glmnet`, are installed in your R environment.
2. Run the script `function.study.name.download_files` by sourcing it in your R console or script environment.
3. The scripts will automatically process the downloaded data and save them in the current working directory as `BRCA_data_meth.csv`, `BRCA_data_meth_range.csv`, and `need_gene_exp_13988(13982).csv`.

### Running Analysis Models

After preparing the data files, you can proceed to run various predictive models to analyze the data further. Here are the steps to perform five-fold cross-validation using different modeling approaches:

1. **DeepMethyGene Model**:

   - Run the `DeepMethyGene.py` script in Python to perform five-fold cross-validation and obtain the model's performance metrics. This script uses deep learning techniques to analyze methylation data.

   ```bash
   python DeepMethyGene.py
   ```

2. **Random Forest Model**:

   - Execute the `Random Forest.py` script to apply the Random Forest algorithm on the data. This model will also use five-fold cross-validation to evaluate its accuracy and robustness.

   ```bash
   python Random Forest.py
   ```

3. **SVM Model**:

   - To apply Support Vector Machine (SVM) modeling, run the `SVM.py` script. This script will conduct five-fold cross-validation to assess the performance of the SVM model on your dataset.

   ```bash
   python SVM.py
   ```

4. **geneEXPLORE Model**:

   - Use the `geneEXPLORE.r` script in R to execute a gene expression analysis with five-fold cross-validation. This script focuses on exploring gene expression data to identify potential patterns or significant features.

   ```bash
   Rscript geneEXPLORE.r
   ```

## Dependencies

- python(version 3.8 or higher recommended)
- R (version 3.6 or higher recommended)
- GenomicRanges: For handling genomic intervals and coordinates.
- glmnet: For modeling data using regularized regression.

## Contact Information

For further assistance or any inquiries related to the scripts and their usage, please do not hesitate to contact:

- email:18853818507@163.com

Feel free to reach out with any questions or if you need support with the implementation or troubleshooting of the provided scripts.

