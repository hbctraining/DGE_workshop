---
title: "Gene-level differential expression analysis with DESeq2"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "April 26, 2017"
---

Approximate time: 60 minutes

## Learning Objectives 

* Understanding the different steps in a differential expression analysis in the context of DESeq2
* Exploring different objects in DESeq2 
* Building results tables for comprison of different sample classes

# Differential expression analysis with DESeq2

The final step in the differential expression analysis workflow is fitting the raw counts to the statistical model and performing the statistical test for differentially expressed genes. In this step we essentially want to determine whether the mean expression levels of two different sample groups are significantly different.

<img src="../img/de_theory.png" width="600">

The [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) was published in 2014, but the package is continually updated and available for use in R through Bioconductor. It builds on good ideas for dispersion estimation and use of Generalized Linear Models from the [DSS]() and [edgeR]() methods. 

Differential expression analysis with DESeq2 involves multiple steps as displayed in the flowchart below in blue. These steps are performed when you run the `DESeq()` function and you do not need to perform them separately. Below we have descriptions of each step, since it is important to understand the inner workings of this tool.

<img src="../img/DESeq2_workflow.png" width="500">

## Step 1: Estimation of size factors

<img src="../img/deseq2_workflow_separate_sf.png" width="200">

DESeq2 will automatically estimate the size factors when performing the differential expression analysis. However, if you have already generated the size factors using `estimateSizeFactors()`, as we did earlier, then DESeq2 will use these values.

## Steps 2 & 3: Estimation of variation/dispersion & Shrinkage of gene-wise dispersion estimates

The next step in the differential expression analysis is the estimation of gene-wise dispersions.

<img src="../img/deseq2_workflow_separate_dis.png" width="200">

To accurately model our sequencing counts, we need to generate accurate estimates of within-group variation (variation between replicates of the same samplegroup) for each gene. With only a few (3-6) replicates per group, the estimates of variation for each gene are often unreliable. Therefore, DESeq2 shares information across genes to generate more accurate estimates of variation based on the expression level of the gene using a method called 'shrinkage'. DESeq2 assumes that genes with similar expression levels have similar dispersion or variation of expression and it used the following steps to generate more accurate measures of dispersion:

1. **Estimate the dispersion for each gene**

	To model the dispersion based on expression level (mean normalized counts of replicates), the dispersion for each gene is estimated using maximum likelihood estimation. In other words, **given the normalized count values of the replicates, the most likely estimate of dispersion is calculated**.

2. **Fit a curve to the the gene estimates given expression strength**

	The idea behind fitting a curve to the data is that different genes will have different scales of biological variability, but, over all genes, there will be a distribution of reasonable estimates of dispersion. 

	This curve is displayed as a red line in the figure below, which plots the estimate for the **expected dispersion value for genes of a given expression strength**. Each black dot is a gene with an associated mean expression level and maximum likelihood estimation (MLE) of the dispersion (Step 1).

	<img src="../img/deseq_dispersion1.png" width="400">

3. **Shrink gene-wise dispersion estimates toward the values predicted by the curve**

	The curve allows for more accurate identification of differentially expressed genes when sample sizes are small, and the strength of the shrinkage for each gene depends on :
	
	- how close gene dispersions are from the curve
	- sample size (more samples = less shrinkage)


	**This shrinkage method is particularly important to reduce false positives in the differential expression analysis.** Genes with extremely low levels of variation are shrunken towards the curve, and the more accurate, higher dispersion values are output for fitting of the model and differential expression testing. 

	Dispersion estimates that are slightly above the curve are also shrunk toward the curve for better dispersion estimation; however, genes with extremely high dispersion values are not shrunken toward the curve due to the likelihood that the gene does not follow the modeling assumptions and has higher variability than others for biological or technical reasons [[1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)]. Shrinking the values toward the curve could result in false positives, so these values are not shrunken. These genes are shown surrounded by blue circles below. 

	<img src="../img/deseq_dispersion2.png" width="600">

	> **NOTE:** This is a good plot to examine to ensure your data is a good fit for the DESeq2 model. You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing expression levels. If you see a cloud or a significantly different shape from what we have above, then you should explore your data some more to see if you have contamination (mitochondrial, etc.) or outlier samples.

## Step 4: Generalized Linear Model fit for each gene

The final step in the DESeq2 workflow is fitting the Negative Binomial model for each gene and performing differential expression testing.

<img src="../img/deseq2_workflow_separate.png" width="200">

As discussed earlier, the count data generated by RNA-seq exhibits overdispersion and the statistical distribution used to model the counts needs to account for this overdispersion. DESeq2 uses a negative binomial distribution to model the RNA-seq counts using the equation below:

 <img src="../img/NB_model_formula.png" width="500">
 
DESq2 will use this formula to create the model for each gene, but what we really want to know is the log2 foldchanges between conditions. The log2 foldchanges can be estimated using the normalized counts with the formula:

 <img src="../img/NB_model_formula_betas.png" width="500">

By fitting the model, DESeq2 will determine the **estimates for the log2 foldchanges and their standard error values for each samplegroup relative to the mean expression of all samples**. Generally for NGS count data, there is a large variance associated with the LFC estimates for genes with low read counts, and these weakly expressed genes would be identified as differentially expressed due solely to this variation. To account for this issue and reduce false positives for lowly expressed genes, DESeq2 shrinks the LFC estimates toward zero when a gene low counts and/or high dispersion values among other things.

### Shrunken log2 foldchanges (LFC)

Similar to the previous shrinkage of dispersion estimates, the shrinkage of LFC estimates uses information from all genes to generate more accurate estimates. Specifically, the distribution of LFC estimates for all genes is used (as a prior) to shrink the LFC estimates of "low information" genes. 

<img src="../img/deseq2_shrunken_lfc.png" width="500">

*Illustration taken from the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).*

For example, in the figure above, the green gene and purple gene have the same mean values for the two sample groups (C57BL/6J and DBA/2J), but the green gene has low variation while the purple gene has high levels of variation. For the green gene with low variation, the unshrunken LFC estimate (vertex of the green solid line) is very similar to the shrunken LFC estimate (vertex of the green dotted line), but the LFC estimates for the purple gene are quite different due to the high dispersion. So even though two genes can have similar normalized count values, they can have differing degrees of LFC shrinkage. Notice the LFC estimates are shrunken toward the prior (black solid line).

>**NOTE:** There are scenarios when it is a good idea to turn off LFC shrinkage which happens by default. (More information about this can be found at [http://seqanswers.com/forums/showthread.php?t=49101](http://seqanswers.com/forums/showthread.php?t=49101))*

### Running DESeq2

To run the differential expression pipeline on the **raw counts** in DESeq2, we must first create a DESeqDataSet object as we did in the 'Count normalization' lesson.

```r
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
```

To run the actual differential expression analysis, we then use a single call to the function DESeq() with the DESeqDataSet object as input. By re-assigning the results of the function back to the same variable name (`dds`), we can continue to fill in the `slots` of our `DESeqDataSet` object.

```r
##Run analysis
dds <- DESeq(dds)
```

This function will print out a message for the various steps it performs: 

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
``` 

<img src="../img/deseq2_workflow_separate.png" width="200">


**Everything from normalization to linear modeling was carried out by the use of a single function!** The results of each step were inserted into the object that you initialized.

![deseq1](../img/deseq_obj2.png)

> **NOTE:** There are individual functions available in DESeq2 that would allow us to carry out each step in the workflow in a step-wise manner, rather than a single call. We demonstrated one example (`estimateSizeFactors()`) when generating size factors to create a normalized matrix. 

### Hypothesis testing using the Wald test

The shrunken LFC estimates for each sample group relative to the mean expression of all groups are stored in the DESeqDataSet object. These estimates represent the **model coefficients**, and these coefficients are calculated for all the levels with respect to the design formula and it is regardless of the comparison of interest. 

The model coefficients can be viewed with `coefficients(dds)` to explore the strength of the effect for each factor group relative to the overall mean for every gene. This is great, but usually we are interested in the LFC estimates relative to other sample groups instead of relative to the mean expression of all groups. For pairwise comparisons between 2 sample groups the difference between the LFC estimates between the groups is calculated. DESeq2 uses the **Wald test** to determine if the difference in the log2 fold changes is significantly from zero (the null hypothesis is that this difference is zero). 

We need to specify to DESeq2 which two groups we want to compare. To do this we can use the `contrasts` argument within the `results()` function to extract the results of the Wald test, i.e. perform the DE analysis between 2 sample groups of interest. 

Contrasts can be provided to DESeq2 a couple of different ways:

1. Automatically DESeq2 will use the base factor level of the condition of interest as the base for statistical testing. 
2. Using the `results()` function, specify the factor and it's levels you would like to compare: `results(dds, contrast=c("sex", "F", "M"))`. The level given last is the base level for the comparison.
3. Instead of giving the factor and levels as a vector, you can create a list using the factor levels given in `resultsNames()`. The level given last is the base level for the comparison. For example, if the output of `resultsNames(dds)` is "sexF", "sexM", then you could write the contrast as follows:
	
	```r
	# DO NOT RUN!
	
	contast_sex <- list("sexF", "sexM")
	
	results(dds, contrast=contrast_sex)
	
	```

> **NOTE:** The Wald test can also be used with **continuous variables**. If the variable of interest provided in the design formula is continuous-valued, then the reported log2 fold change is per unit of change of that variable.


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

*Some materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*

***

 

