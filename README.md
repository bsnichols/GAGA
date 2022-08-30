# GAGA: a pipeline for GEM And GWAS Automation

Welcome to the GEM And GWAS Automisation (GAGA) Pipeline!

The GAGA pipeline is an R package for GWAS And GEM Automisation. The pipeline has been designed to provide a user-friendly method of analysing multiple traits using both GEM and GWAS analysis automatically in hours. GAGA relies on the pre-existing GAPIT3 package[^Wang] to complete the GWAS analysis. 

In this repo you will find the  GAGA pipeline script, demo data and the documentation in this README file.

![alt text](https://github.com/bsnichols/GAGA/blob/main/paper/pipeline_flowchart.png "Flowchart")

The pipeline can analyse multiple traits from a single file and output the results for each trait in individual folders, with a collated results summary for ease of interpretation. The pipeline consists of two main stages. In the first, named Data management, the multiple traits are separated and organised into unique folders. If the inputted data requires transformation (as some statistical analyses insist on normality), this can optionally be included to compare the raw, log-transformed and square-rooted data using ANOVA. The most normally distributed data is processed using a linear mixed model and the estimated means are saved for the next stage. The second stage consists of the main Analyses. GEM and/or GWAS analysis can be performed. During the GWAS analysis, it is possible to select the types of models. The pipeline uses Least Sum of Squares (LSS) to determine the best-fit model out of those selected and carries this forward into plotting. Once the analysis/analyses are completed, Q values and a False Discovery Rate (FDR) are calculated. This information, along with the results are plotted using the Manhattan plotter and saved. The colour of the Manhattan plots is customisable. This stage is repeated for all the traits from the multiple traits file. Once the pipeline has carried out the Analyses stage on all the traits, the pipeline will summarise the results by creating a table of traits with significant GEM and/or GWAS results.

While the GEM analysis and plotter functions of the pipeline are currently configured for the A/C genome requirements of *Brassica napus*, this can be readily adapted for other species.

## 1 Setting up

The pipeline runs within the folder it is stored in. To begin, place GAGA.R in the folder containing the data required by the pipeline. This should include:

- the trait input data in the .csv format, with the prefix 'totest_', containing the traits to be run (more information on this is detailed in 2)
- data for GWAS (Fell *et al.*, preprint[^Fell])
    * q matrix in the .txt format, with the prefix 'qmatrix_'
    * SNP file in the .hmp.txt format, with the prefix 'snpfile_'
- data for GEM (from Woolfenden, 2022[^Woolfenden])
    * rpkm_AT2018 in the .txt format, with the prefix 'unigenefile_'
    * Marker to At_AT2018 in the .csv format, with the prefix 'marker_'
- and data for plotting (from Woolfenden, 2022[^Woolfenden])
    * DirectionsAC_AT2018 in the .tsv format, with the prefix 'directions_'
    * Ath_Mapping in the .tsv format, with the prefix 'mapping_'

Take a look at the demo data in this repo for how to name and format your data.

When you open up GAGA.R, there are six variables you will need to manually set before running the pipeline. ```setwd()``` on line 10 requires you to put the path to the directory containing GAGA.R and the data files listed above. 

```runstats``` needs to be set to ```TRUE``` or ```FALSE``` depending on whether or not you wish to transform your data and ```colno``` is dependent on the type of Linear Mixed Model you wish to run, should you be transformin your data. For a further explanation on ```runstats``` and ```colno```, refer to 1.1.

```rungem``` runs the GEM analysis and ```rungwas``` runs the GWAS analysis and these need to be set to ```TRUE``` or ```FALSE``` depending on which you wish to run on your data. If you are running a GWAS, the ```gwasmodels``` should be adapted for the models you wish to run in GAPIT3[^Wang].

To run the pipeline, ensure these variables are set and that the correct data files in the directory, and then the pipeline will automatically run when sourced.

### 1.1 Installing GAPIT package

This pipeline relies on GAPIT3. Documentation for how to install GAPIT3 can be found on the GAPIT website: https://zzlab.net/GAPIT/.

I have found the best way to download it and its dependencies to be the following:

```R
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3")
library(GAPIT3)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

library(BiocManager)
BiocManager::install("multtest")
BiocManager::install("scatterplot3d")
BiocManager::install("gplots")
```

## 2 Data management

The input data format is dependent on whether or not the data has been transformed. 

For already transformed data, the format is as follows:

| genotype_id | traitone | traittwo | traitthree | trait... |
|-|-|-|-|-|

For data to be transformed, the format is as follows:

| genotype_id | location | rep | traitone | traittwo | traitthree | trait... |
|-|-|-|-|-|-|-|

See the data files in the traitdataformats folder in this repository for help.

We recommend that you use lowercase letters and no numbers or symbols for your trait names.

As standard, the pipeline for data to be transformed it assumes that genotype_id, location and rep feature before the traits. If your data is in a different format to this, you will need to adjust the transformation stage to match using the guide in 2.1.

### 2.1 Transforming the data

To transform the data, set ```runstats``` on line 11 to ```TRUE```. 

The pipeline is currently set up to run the following as a Linear Mixed Model (LMM) on line 154 __SUBJECT TO CHANGE__: 

```R
LMMmod<-lmer(trait~(1|location)*genotype_id, data = MyDataframe)
```

If your data requires a different model, replace ```trait~(1|location)*genotype_id``` with the model of your choice. If the number of columns in your input data is therefore different from the format in 2, replace ```colno``` on line 12 with the number of columns that come before your first trait. For example, in the format shown there are 3 columns (genotype_id, location and rep) before the first trait and so ```colno = 3```.

## 3 Manhattan plots

In the Manhattan plots, the significant Q-value threshold line is drawn in blue and the significant False Discovery Rate (FDR) is drawn in red.

The Manhattan plotter features twice within the code; once within the GEM analysis and once within the GWAS analysis. To change the colours of points in the Manhattan plots, you will need to alter ```Chrom_colors``` to the colours of your choice in both regions within the script. As standard, the Manhattan plots produce a rainbow Manhattan plot of multiple colours. 

[^Fell]: Fell, H., *et al.* (2022) "Novel gene loci associated with susceptibility or cryptic quantitative resistance to Pyrenopeziza brassicae in Brassica napus.", *preprint*.
[^Wang]: Wang J., Zhang Z. (2021) ‘GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction, Genomics, Proteomics & Bioinformatics’, doi: https://doi.org/10.1016/j.gpb.2021.08.005.
[^Woolfenden]: Woolfenden, H. (2022) ‘Pyrenopeziz Resistance project’ Github repository, doi: https://doi.org/10.5281/zenodo.6546233.
