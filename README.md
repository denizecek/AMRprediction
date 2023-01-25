 
# AMRprediction

This repository contains scripts for predicting AMR class (R-S) using nucleotide k-mers, amino acid k-mers and SNP features. The overiew of each script (in order) is presented below, details are provided as comments within the R script files.

## Feature Processing

The [feature processing script](Feature_Processing.R) are used to:

- Generate nucleotide k-mer counts from contig.fasta files
- Extract SNP information from VCF files
- Generate AA k-mer counts from protein.fasta files

## Feature Selection

The [feature selection script](Feature_Selection.R) are used to:

- Perform feature elimination (keep samples with absolute mean difference less than a set threshold between the 3 AMR classes)
- Feature selection using [Boruta](https://doi.org/10.18637/jss.v036.i11)

## Machine Learning Models

The [machine learning model scripts](Machine_Learning_Models.R) are used to:

- Split the data into train (80%), validation (10%) and test (10%) datasets
- Train and evaluate:
	- Random Forest models
	- Support Vector Machine (SVM) models
	- Stochastic Gradient Boosting (GBM) models
	- eXtreme Gradient Boosting (xGBoost) models
