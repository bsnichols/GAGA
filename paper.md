---
title: 'GAGA: An R package for high-speed GWAS And GEM Automisation'
tags:
  - R
  - GWAS
  - GEM
  - Phenotyping
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
- name: Rachel Wells
    orcid: 0000-0002-1280-7472
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
affiliations:
 - name: John Innes Centre, UK
   index: 1
 - name: University College London, UK
   index: 2
date: 11 July 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Over the last decade, high-throughput phenotyping by drones and phenotyping platforms have improved phenotyping capabilities in the field and lab. The accuracy of high-throughput phenotyping methods is also improving, with new trait extraction techniques being developed rapidly for structural levels whole canopy to microscopic. With this influx of trait extraction has come a backlog of genomic association analysis. Traditionally, single traits are compared to genomic information such as gene expression or single nucleotide polymorphisms using Genome Wide Association Study (GWAS) analysis or Gene Expression Marker (GEM) analysis. As it can take several days to analyse a single trait using these methods, the phenotyping bottleneck has progressed from trait extraction from images to genomic association. While the tools to perform these association studies are available in many forms, the manual process of implementing them is causing the backlog.  

# Statement of need

The `GAGA` pipeline is an R package for GWAS And GEM Automisation. The pipeline has been designed to provide a user-friendly method of analysing multiple traits using both GWAS and GEM analysis automatically in hours. GAGA relies on the pre-existing [@GAPIT3] package to complete the GWAS analysis [@Wang:2020]. The design of the pipeline can be broken down into three stages: transformation, GEM analysis and GWAS analysis \autoref{fig:flowchart}. The pipeline is designed to take in either raw replicate data and transform it for analysis, or ready transformed data. The pipeline then performs a GEM analysis on the data and outputs a Manhattan plot. It then performs a GWAS using GAPIT and can be adjusted to run any of the models available in the GAPIT package. The best-fit model feature compares the results from the models run in the GWAS and outputs a Manhattan plot for the model that best fits the data. Outputs from the pipeline are organised into individual trait folders and an overview of all traits analysed is available to help quickly navigate the results. 

Each of the steps performed by GAGA are currently run manually and can take several hours per trait. GAGA was designed to automate the genomic association analysis pipeline being manually performed by the [@BRAVO] (Brassica, Rapeseed And Vegetable Optimisation) project, and can run multiple traits automatically, with each trait taking between 20-40 minutes, dependent on the speed of your machine. The development of the pipeline has lead to the discovery of a potential pathogen related gene in *Brassica napus* [@Jacott:2022]. The pipeline has also been used in other phenotyping studies related to the [@BRAVO] project in seed and pod anlysis, branching patterns, flowering time and whole plant development at the John Innes Centre, the University of York, the University of Aberystwyth and Rothamsted Research. While the pipeline has thus far been used to analyse *Brassice napus* only, it can be tailored to study other species. The source code for GAGA with documentation and demo data has been archived to Zenodo with the linked DOI: [LINK TO DOI].

# Figures

![Flowchart for the GAGA pipeline.\label{fig:flowchart}](figure.png)

# Acknowledgements

BRAVO funding, Turing and BBSRC funding
Catherine Jacott for her contributions during the testing of this software.

# References

