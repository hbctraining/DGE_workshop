**THIS REPO IS ARCHIVED, PLEASE GO TO [https://hbctraining.github.io/main](https://hbctraining.github.io/main) FOR CURRENT LESSONS.**

# Differential gene expression workshop

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R/) | 1.5-day workshop (~10 hours of trainer-led time)|

### Description

This repository has teaching materials for a **1.5-day**, hands-on **Introduction to differential gene expression (DGE) analysis** workshop. The workshop will lead participants through performing a differential gene expression analysis workflow on RNA-seq count data using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/).

### Learning Objectives

- QC on count data using Principal Component Analysis (PCA) and hierarchical clustering
- Using DESeq2 to obtain a list of significantly different genes
- Visualizing expression patterns of differentially expressed genes
- Performing functional analysis on gene lists with R-based tools

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.

### Lessons

Below are links to the lessons and suggested schedules:

* [Click here for schedule using Salmon count matrix](https://hbctraining.github.io/DGE_workshop_salmon/schedule/)
* [Click here for schedule using FeatureCounts count matrix](schedule/1.5-day.md)


### Installation Requirements

1. Download the most recent versions of R and RStudio for your laptop:

 - [R](https://cran.r-project.org/) 
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

2. Install the following packages using the instructions provided below.

> **NOTE:**  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable but know that it can take awhile.

(a) Install the below **packages** on your laptop **from CRAN**. You DO NOT have to go to the CRAN webpage; you can use the following function to install them one by one:

```r
install.packages("insert_first_package_name_in_quotations")
install.packages("insert__second_package_name_in_quotations")
& so on ...
```

Packages to install from CRAN (note that these package names are case sensitive!):

* BiocManager
* RColorBrewer
* pheatmap
* ggrepel
* devtools
* tidyverse


(b) Install the below **packages from Bioconductor**, using `BiocManager::install()` function 7 times for the 7 packages:

```
BiocManager::install("insert_first_package_name_in_quotations")
BiocManager::install("insert_second_package_name_in_quotations") 

```

Packages to install from Bioconductor (note that these package names are case sensitive!):

* DESeq2
* clusterProfiler
* DOSE
* org.Hs.eg.db
* pathview
* DEGreport
* EnsDb.Hsapiens.v86
* AnnotationHub
* ensembldb


3. Finally, please check that all the packages were installed successfully by loading them one at a time using the library() function.  

```r
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
```

4. Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```



****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
