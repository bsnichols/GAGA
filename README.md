# The GAGA Pipeline

Welcome to the GEM And GWAS Automisation (GAGA) Pipeline!

The GAGA pipeline is an R package for GWAS And GEM Automisation. The pipeline has been designed to provide a user-friendly method of analysing multiple traits using both GEM and GWAS analysis automatically in hours. GAGA relies on the pre-existing GAPIT3 package to complete the GWAS analysis, and the GEM script developed by Woolfenden (2020). 

In this repo you will find the  GAGA pipeline script, demo data and a tutorial in the README below.

The design of the pipeline can be broken down into three stages: data organisation, GEM analysis and GWAS analysis section. The pipeline is designed to take in either raw replicate data and transform it for analysis, or ready transformed data. The pipeline then performs a GEM analysis on the data and outputs a Manhattan plot. It then performs a GWAS using GAPIT and can be adjusted to run any of the models available in the GAPIT package. The best-fit model feature compares the results from the models run in the GWAS and outputs a Manhattan plot for the model that best fits the data. Outputs from the pipeline are organised into individual trait folders and an overview of all traits analysed is available to help quickly navigate the results.

## 1 Data management

The pipeline runs within the folder it is stored in. To begin, place GAGA.R in the folder containing the data required by the pipeline. This should include:

- the trait input data, with the prefix 'totest_', containing the traits to be run (more information on this is detailed in 1.1)
- data for GEM
    * rpkm_AT2018 in the .txt format
    * Marker to At_AT2018 in the .csv format
- data for GWAS
    * q matrix in the .txt format, with the prefix 'qmatrix_'
    * SNP file in the .hmp.txt format, with the prefix 'snpfile_'
- and data for plotting
    * DirectionsAC_AT2018 in the .tsv format
    * Ath_Mapping in the .tsv format

__WHAT ARE THE OTHER FILES CALLED__

Take a look at the demo data in this documentation for how to name and format your data.

When you open up GAGA.R, there are three variables you will need to set. ```setwd()``` on line 10 requires you to put the path to the directory containing GAGA.R and the data files listed above. ```runstats``` needs to be set to ```TRUE``` or ```FALSE``` depending on whether or not you wish to transform your data and ```colno``` is dependent on the type of Linear Mixed Model you wish to run, should you be transformin your data. For a further explanation on ```runstats``` and ```colno```, refer to 1.1.

### 1.1 Trait input data

The input data format is dependent on whether or not the data has been transformed. 

For already transformed data, the format is as follows:

| genotype_id | traitone | traittwo | traitthree | traitn |
| -------------|:-------------:|:-------------:|:-------------:| -----:|
| | | | | | |

For data to be transformed, the format is as follows:

| genotype_id | location | rep | traitone | traittwo | traitthree | traitn |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|:-------------:| -----:|
| | | | | | |

As standard, the pipeline for data to be transformed it assumes that genotype_id, location and rep feature before the traits. If your data is in a different format to this, you will need to adjust the transformation stage to match using the guide in 1.2.

### 1.2 Transforming the data

To transform the data, set ```runstats``` on line 11 to ```TRUE```. 

The pipeline is currently set up to run the following as a Linear Mixed Model (LMM) on line 154 __SUBJECT TO CHANGE__: 

```R
LMMmod<-lmer(trait~(1|location)*genotype_id, data = MyDataframe)
```

If your data requires a different model, replace ```trait~(1|location)*genotype_id``` with the model of your choice. If the number of columns in your input data is therefore different from the format in 1.1, replace ```colno``` on line 12 with the number of columns that come before your first trait. For example, in the standard format there are 3 columns before the first trait and so ```colno = 3```.

