---
title: "Gene-level differential expression analysis with DESeq2"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "April 26, 2017"
---

Approximate time: 60 minutes

## Learning Objectives 

* 1
* 2
* 3



## Differential expression analysis with DESeq2

The final step in the differential expression analysis is the actual fitting of the raw counts to the statistical model and testing for differentially expressed genes. Essentially we want to determine whether the mean expression levels of two different sample groups are significantly different.

<img src="../img/de_theory.png" width="600">

DESeq2 builds on good ideas for dispersion estimation and use of Generalized Linear Models from the DSS and edgeR methods. The [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) was published in 2014, but the package is continually updated and available for use in R through Bioconductor.

Differential expression analysis with DESeq2 requires multiple steps, as displayed below.

<img src="../img/DESeq2_workflow.png" width="500">

### Design formula

When performing differential expression analysis, it is a good idea to know what sources of variation are present in your data, either by exploration during the QC or prior knowledge. Once you know the major sources of variation, you can remove them prior to analysis or control for them in the statistical model. For example, if you know that animal sex or age is a significant source of variation in your data, then it needs to be included in your model. The **design formula** or **model formula** needs to have **all of the factors in your metadata that account for major sources of variation** in your data. The last factor entered in the formula should be the condition of interest. 

For example, suppose you have the following metadata appears as follows:

If you want to examine the expression differences between treatments, and you know that major sources of variation include 'sex' and 'age', then your design formula would be:

`design <- ~ sex + age + treatment`

***
**Exercises**

1. Suppose you wanted to study the expression difference between the two age groups, and major sources of variation were 'sex' and 'treatment', how would the design formula be written?
2. Based on our `metadata` dataframe, which factors could we include in our design formula?
3. What would you do if you wanted to include a factor in your design formula that is not in your metadata? 

***


#### Complex designs

DESeq2 also allows for the analysis of complex designs. You can explore interactions or difference of differences by specifying for it in the design formula. For example, if you wanted to explore the effect of sex on the treatment  effect, you could specify for it in the design formula as follows: 

`design <- ~ sex + age + treatment + sex:treatment`

Since the interaction term `sex:treatment` is last in the formula, the results output from DESeq2 will output results for this term. Alternatively, as recommended in the [DESeq2 vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf), we could create a new factor variable in our metadata based on the two interaction factors (i.e. "Ftreated", "Fcontrol", "Mtreated", "Mcontrol"). 


### Estimate size factors

After you have your design formula saved to a variable, you can input the **design formula**, in addition to the **raw counts** and **metadata** into the `DESeq()` function to perform the steps in the DESeq2 differential expression analysis. **You must input the RAW counts (not normalized) for the analysis.**

<img src="../img/deseq2_workflow_separate_sf.png" width="200">

The first step in the differential expression is to estimate the size factors, which is exactly what we already did to normalize the raw counts. DESeq2 will automatically estimate the size factors when performing the differential expression analysis if you haven't already done so. However,  you have already generated the size factors, then DESeq2 will use these values. 

### Estimate variation / dispersion

The next step in the differential expression analysis is the estimation of gene-wise dispersions.

<img src="../img/deseq2_workflow_separate_dis.png" width="200">

To accurately model our sequencing counts, we need to generate accurate estimates of within-group variation (variation between replicates of the same samplegroup) for each gene. With only a few (3-6) replicates per group, the estimates of variation for each gene are often unreliable. Therefore, DESeq2 shares information across genes to generate more accurate estimates of variation based on the expression level of the gene using a method called 'shrinkage'. DESeq2 assumes that genes with similar expression levels have similar dispersion or variation of expression. DESeq2 generates more accurate measures of dispersion using the following steps:

1. **Estimate the dispersion for each gene separately**

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

	> **NOTE:** This is a good plot to examine to ensure your data is a good fit for the DESeq2 model. You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing expression levels. If you see a cloud or different shapes, then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.

### Generalized Linear Model fit for each gene

The final step in the DESeq2 workflow is fitting the Negative Binomial model for each gene and performing differential expression testing.

<img src="../img/deseq2_workflow_separate.png" width="200">

As discussed earlier, the count data generated by RNA-Seq exhibits overdispersion and the statistical distribution used to model the counts needs to account for this overdispersion. DESeq2 uses a negative binomial distribution to model the RNA-Seq counts using the equation below:

 <img src="../img/NB_model_formula.png" width="500">
 
DESq2 will use this formula to create the model for each gene, but what we really want to know is the log2 foldchanges between conditions. The log2 foldchanges can be estimated using the normalized counts with the formula:

 <img src="../img/NB_model_formula_betas.png" width="500">

By fitting the model, DESeq2 will determine the **estimates for the log2 foldchanges and their standard error values for each samplegroup relative to the mean expression of all samples**. 

#### Shrunken log2 foldchanges (LFC)

Generally for NGS count data, there is a large variance associated with the LFC estimates for genes with low read counts, and these weakly expressed genes would be identified as differentially expressed due solely to this variation. To account for this issue and reduce false positives for lowly expressed genes, DESeq2 shrinks the LFC estimates toward zero when the information for a gene is low, which could include:

- Low counts
- High dispersion values

Similar to the previous shrinkage of dispersion estimates, the shrinkage of LFC estimates uses information from all genes to generate more accurate estimates. Specifically, the distribution of LFC estimates for all genes is used (as a prior) to shrink the LFC estimates of genes with little information or high dispersion toward more likely (lower) LFC estimates. 

<img src="../img/deseq2_shrunken_lfc.png" width="500">

*Illustration taken from the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).*

For example, in the figure above, the green gene and purple gene have the same mean values for the two sample groups (C57BL/6J and DBA/2J), but the green gene has low variation while the purple gene has high levels of variation. For the green gene with low variation, the unshrunken LFC estimate (vertex of the green solid line) is very similar to the shrunken LFC estimate (vertex of the green dotted line), but the LFC estimates for the purple gene are quite different due to the high dispersion. So even though two genes can have similar normalized count values, they can have differing degrees of LFC shrinkage. Notice the LFC estimates are shrunken toward the prior (black solid line).


>**NOTE:** If very large expected fold changes for a number of individual genes are expected, but not so many large fold changes that the width of the prior adjusts to allow such large fold changes, then you may want to turn off LFC shrinkage.
>
>The reason is that shrinking of fold changes requires that the software can estimate the range of reasonable values for LFC by looking at the distribution of LFCs (particularly the upper quantile of the absolute LFC). But there might be precise fold changes which are above this upper quantile, and so the prior is too narrow for the targeted genes. The prior might then be a bad assumption for this type of dataset, so it's reasonable to turn it off. *(Response from Mike Love, creator of DESeq2 (http://seqanswers.com/forums/showthread.php?t=49101))*


#### Hypothesis testing using the Wald test

The shrunken LFC estimates are given for each sample group relative to the mean expression of all groups. These estimates represent the **model coefficients**, and these coefficients are calculated regardless of the comparison of interest. The model coefficients can be viewed with `coefficients(dds)` to explore the strength of the effect for each factor group relative the overall mean for every gene. 

Generally we are interested in the LFC estimates relative to other sample groups instead of to the mean expression of all groups. To do this, we must test if the difference in the log2 fold changes between groups is zero. To determine whether the difference in shrunken LFC estimates differs significantly from zero, the **Wald test** is used. The Wald test is generally used to make pair-wise comparisons (i.e. compare the LFCs from two different conditions).

To indicate to DESeq2 the two groups we want to compare, we can use **contrasts** to perform differential expression testing using the Wald test. 

Contrasts can be provided to DESeq2 a couple of different ways:

1. Automatically DESeq2 will use the base factor level of the condition of interest as the base for statistical testing. 
2. Using the `results()` function, specify the factor and it's levels you would like to compare: `results(dds, contrast=c("sex", "F", "M"))`. The level given last is the base level for the comparison.
3. Instead of giving the factor and levels as a vector, you can create a list using the factor levels given in `resultsNames()`. The level given last is the base level for the comparison. For example, if the output of `resultsNames(dds)` is "sexF", "sexM", then you could write the contrast as follows:
	
	```r
	
	# DO NOT RUN!
	
	contast_sex <- list("sexF", "sexM")
	
	results(dds, contrast=contrast_sex)
	
	```

#### Multiple test correction

If we used the p-value directly from the Wald test with a significance cut-off of 0.05 (α = 0.05), then 5% of all genes would be called as differentially expressed (i.e. 5% False positive genes). The more genes we test, the more 'false positives' we discover. DESeq2 helps reduce the number of genes tested by removing those genes unlikely to be significantly DE, such as those with low number of counts and outlier samples. However, we still need to correct for multiple testing, and there are a few common approaches:

- **Bonferroni:** Reject any hypothesis with p-value ≤ α/m. **This is a very conservative approach with a high probability of false negatives.**
- **FDR / Benjamini-Hochberg:** Rank j / m multiplied by the FDR levels. This approach is designed to control the proportion of false positives among the set of rejected hypotheses
- **Q-value:** The minimum FDR that can be attained when calling that feature significant. For example, if gene X has a q-value of 0.013 it means that 1.3% of genes that show p-values at least as small as gene X are false positives

In DESeq2, the p-values attained by the Wald test are corrected for multiple testing using the Benjamin and Hochberg method. The p-adjusted values should be used to determine significant genes.

 




***
# Turning off prior
 You can turn off the LFC shrinkage using `DESeq(dds, betaPrior=FALSE)`
 
# Model matrix
In DESeq2 and most other DE tools, you will assign your samples to specific conditions using a 'model matrix' or 'design matrix'. For example, we can assign our samples to conditions in the model matrix using binary (0,1) notation:

|  | ctrl | oe | kd |
| ----- |:-----:|:-----:|:-----:|
| sample1 | 1 | 0 | 0 |
| sample2 | 1 | 0 | 0 |
| sample3 | 0 | 1 | 0 |
| sample4 | 0 | 1 | 0 |
| sample5 | 0 | 0 | 1 |
| sample6 | 0 | 0 | 1 |

Samples 1 and 2 are controls, samples 3 and 4 are overexpression, and samples 5 and 6 are knockdown. This information is utilized to inform the model about which replicates should be used to estimate the **log2 foldchanges (LFC)** for each level of each factor in the model.

>For example, sample1 and sample2 should be used to estimate the LFC for the expression of the control group relative to mean expression of all groups, and the model formula for this comparison would be:
>
><img src="../img/NB_model_formula_betas_example3.png" width="400">
