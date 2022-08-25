# The GAGA Pipeline

Welcome to the GEM And GWAS Automisation (GAGA) Pipeline!

The GAGA pipeline is an R package for GWAS And GEM Automisation. The pipeline has been designed to provide a user-friendly method of analysing multiple traits using both GEM and GWAS analysis automatically in hours. GAGA relies on the pre-existing GAPIT3 package to complete the GWAS analysis, and the GEM script developed by Woolfenden (2020). 

In this repo you will find the  GAGA pipeline script, demo data and a tutorial in the README below.

The design of the pipeline can be broken down into three stages: data organisation, GEM analysis and GWAS analysis section. The pipeline is designed to take in either raw replicate data and transform it for analysis, or ready transformed data. The pipeline then performs a GEM analysis on the data and outputs a Manhattan plot. It then performs a GWAS using GAPIT and can be adjusted to run any of the models available in the GAPIT package. The best-fit model feature compares the results from the models run in the GWAS and outputs a Manhattan plot for the model that best fits the data. Outputs from the pipeline are organised into individual trait folders and an overview of all traits analysed is available to help quickly navigate the results.

## 1 Data management

### 1.1 Data input

The input data format is dependent on whether or not the data has been transformed. 

For already transformed data, the format is as follows:

| genotype_id | traitone | traittwo | traitthree | trait... |
| ------------- |:-------------:|:-------------:|:-------------:| -----:|
| | | | |

For data to be transformed, the format is as follows:

| genotype_id | location | rep | traitone | traittwo | traitthree | trait... |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|:-------------:| -----:|
| | | | | | |

As standard, the pipeline for data to be transformed it assumes that genotype_id, location and rep feature before the traits. If your data is in a different format to this, you will need to adjust the transformation stage to match using the guide in 1.2.

### 1.2.1 Transforming the data

The pipeline is currently set up to run the following as a Linear Mixed Model (LMM) on line 154 __SUBJECT TO CHANGE__: 

```R
LMMmod<-lmer(trait~(1|location)*genotype_id, data = MyDataframe)
```
.

If your data requires a different model, replace ```R trait~(1|location)*genotype_id ``` with the model of your choice. If the number of columns in your input data is therefore different, replace ```R colno ``` with the number of columns that comes before your first trait.

#### 1.2.1 Adjusting the transformation stage

