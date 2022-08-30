---
title: 'GAGA: A pipeline for GEM and GWAS Automation'
tags:
  - R
  - GWAS
  - GEM
  - Phenomics
  - Genomics
authors:
  - name: Bethany S. Nichols
    orcid: 0000-0003-1280-2070
    corresponding: true
    affiliation: 1
  - name: Hugh Woolfenden
    orcid: 0000-0002-9800-7862
    equal-contrib: true 
    affiliation: 1
- name: Jo Hepworth
    orcid: 0000-0002-4621-8414
    equal-contrib: true 
    affiliation: 1
- name: Joshua Williams
    orcid: 0000-0002-8423-5549
    equal-contrib: true 
    affiliation: 2
- name: Lars Østergard
    orcid: 0000-0002-8497-7657
    equal-contrib: true 
    affiliation: 1
- name: Richard Morris
    orcid: 0000-0003-3080-2613
    equal-contrib: true 
    affiliation: 1
 - name: Rachel Wells
    orcid: 0000-0002-1280-7472
    equal-contrib: true 
    affiliation: 1
affiliations:
 - name: John Innes Centre, UK
   index: 1
 - name: University College London, UK
   index: 2
date: 26 August 2022
bibliography: paper.bib

# Summary

Over the last decade, high-throughput phenotyping has transformed capabilities in the field and lab [@Mir:2019]. The accuracy of high-throughput phenotyping methods has improved steadily, with new trait extraction techniques of increasing resolution being developed for spatial scales ranging from microscopic to whole canopy and field level. With this surge in trait extraction has come a backlog of genomic association analysis [@Nagy:2020]. Traditionally, single traits are compared to genomic information such as gene expression or single nucleotide polymorphisms using Gene Expression Marker (GEM) analysis or Genome Wide Association Study (GWAS) analysis [@Chen:2014]. As it can take several days to analyse a single trait using these methods, the phenotyping bottleneck has progressed from trait extraction from images to genomic association. While the tools to perform these association studies are readily available, the manual processing and decision-making required to carry out all associated steps can be time-consuming and relies on expert domain knowledge. Here, we describe an automated pipeline that carries out all associated steps from data processing to trait association and analysis without the need for human intervention.


# Statement of need

We describe the `GAGA` pipeline for **G**EM **A**nd **G**WAS **A**utomation \autoref{fig:flowchart}. The pipeline is written in RStudio [@RTeam:2020] and has been designed to provide a fast, user-friendly method for analysing multiple traits using both GEM analysis and GWAS. GAGA requires the installation of `GAPIT3` for GWAS analysis [@Wang:2020]. The design of the pipeline can be broken down into two stages: Data management and Analyses. The pipeline takes in either raw replicate data and transform it for analysis, or ready transformed data. As it can process multiple traits, during Data management the pipeline separates the traits and organises them into unique folders. The pipeline then performs selected Analyses on the traits. Depending on user selections (and available data) it will run a GEM analysis and/or a GWAS analysis. GEM analyses are carried out using an adapted script [@Fell:2022]. Following a GEM analysis, the pipeline outputs a Manhattan plot labelled with significance markers. If a GWAS analysis is selected, the pipeline can be configured to run any of the models available in the GAPIT package. The best-fit model feature compares the results from the models run in the GWAS and outputs a Manhattan plot, labelled with significance markers, for the model that best fits the data. Outputs from the pipeline are organised into individual trait folders and an overview of all traits analysed is available to help quickly navigate the results. 

GAGA implements practical domain knowledge to mimic the steps that would otherwise need to be performed manually. Running such an analysis manually requires an experienced researcher and can take several hours per trait. GAGA was developed to automate these steps for multiple traits. In the current version each trait takes between 20-40 minutes, depending on available processing power. The pipeline has been used in several *Brassica napus* phenotyping studies, including seed and pod analysis, branching patterns, flowering time and whole plant development at the John Innes Centre, the University of Nottingham, the University of Aberystwyth and Rothamsted Research. As the pipeline can process large numbers of traits automatically, this has enabled the rapid comparison of multiple disease association studies in *B. napus* [@Jacott:2022], allowing an association of a region controlling disease resistance to be discovered which would likely have remained undetected using manual analysis methods. While the GEM aspect of the pipeline is currently configured for the A/C genome requirements of *B. napus*, this can be readily adapted for any species. 

The source code for GAGA with documentation and demo data has been archived to Zenodo with the linked DOI: [LINK TO DOI]. 

# Figures

![A pipeline for GEM And GWAS Automation (GAGA). The pipeline can analyse multiple traits from a single file and output the results for each trait in individual folders, with a collated results summary for ease of interpretation. The pipeline consists of two main stages. In the first, named Data management, the multiple traits are separated and organised into unique folders. If the inputted data requires transformation (as some statistical analyses insist on normality), this can optionally be included to compare the raw, log-transformed and square-rooted data using ANOVA. The most normally distributed data is processed using a linear mixed model and the estimated means are saved for the next stage. The second stage consists of the main Analyses. GEM and/or GWAS analysis can be performed. During the GWAS analysis, it is possible to select the types of models. The pipeline uses Least Sum of Squares (LSS) to determine the best-fit model out of those selected and carries this forward into plotting. Once the analysis/analyses are completed, Q values and a False Discovery Rate (FDR) are calculated. This information, along with the results are plotted using the Manhattan plotter and saved. The colour of the Manhattan plots is customisable. This stage is repeated for all the traits from the multiple traits file. Once the pipeline has carried out the Analyses stage on all the traits, the pipeline will summarise the results by creating a table of traits with significant GEM and/or GWAS results.\label{fig:flowchart}](https://github.com/bsnichols/GAGA/blob/main/pipeline_flowchart.png)

# Acknowledgements

We thank our colleague Catherine Jacott from the John Innes Centre for her contributions during the testing of this software. LO and RW acknowledge funding from the BBSRC Institute Strategic Programme 'Genes in the environment’ (BB/P013511/1), RM, JH and RW acknowledge support from the BBSRC strategic LoLa ‘Brassica Rapeseed and Vegetable Optimisation’ (BB/P003095/1), RW, LO, RM and BN acknowledge support from the BBSRC's Catalyst Partnership in Artificial Intelligence between the Alan Turing Institute and the Norwich Bioscience Institute (BB/V509267/1). 

# References

