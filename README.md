---
Audience: Biologists
Computational Skills: Intermediate R
Prerequisites: Introduction to R
Duration: 1-day workshop (~6.5 hours of trainer-led time)
---


# Differential gene expression workshop

This repository has teaching materials for a **1-day**, hands-on **Introduction to differential gene expression (DGE) analysis** workshop. The workshop will lead participants through performing a differential gene expression analysis workflow on RNA-seq count data using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://github.com/hbctraining/Intro-to-R).

### Learning Objectives

- QC on count data using Principal Component Analysis (PCA) and heirarchical clustering
- Using DESeq2 to obtain a list of significantly different genes
- Visualizing expression patterns of differentially expressed genes
- Performing functional analysis on gene lists with R-based tools

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.

### Contents

#### Differential Gene Expression (DGE) using RNA-seq raw counts data
| Lessons            | Duration |
|:------------------------|:----------|
|[Setting up and DGE overview](lessons/10_DGE_setup_and_overview.md) | 70 min |
|[Introduction to count normalization](lessons/11_DGE_count_normalization.md) | 60 min |
|[QC using principal component analysis (PCA) and heirarchical clustering](lessons/12_DGE_QC_analysis.md) | 90 min |
|[Getting started with DESeq2](lessons/13_DGE_DESeq2_analysis.md) | 70 min |
|[Pairwise comparisons with DEseq2](lessons/14_DGE_DESeq2_analysis2.md) | 45 min |
|[Visualization of DGE analysis results](lessons/15_DGE_visualizing_results.md) | 45 min |
|[Complex designs with DESeq2 (LRT)](lessons/16_DGE_LRT.md) | 60 min |
|[BiomaRt and Ensembl](lessons/17_Ensembl_biomart.md) | 30 min |
|[Functional Analysis](lessons/18_functional_analysis.md) | 85 min |

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
