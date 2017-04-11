# Differential gene expression analysis with DESeq2

The goal of RNA-Seq is often to perform differential expression testing to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the condition(s) of interest. 

To determine the expression levels of genes, an RNA-Seq workflow is followed with the steps detailed in the image below. All steps are performed on the command line (Linux/Unix) through the generation of the read counts per gene. The differential expression analysis and any downstream functional analysis are generally performed in R using R packages specifically designed for the complex statistical analyses required to determine whether genes are differentially expressed.

<img src="../img/rnaseq_full_workflow.png" width="400">


The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher level of expression of that gene in the sample. 

<img src="../img/deseq_counts_overview.png" width="500">




**The goal of differential expression analysis is to determine whether the variation in expression (counts) between groups is significantly greater than variation within groups (replicates).** To test for significance, we need an appropriate statistical model that accurately performs normalization (to account for differences in sequencing depth, etc.) and variance modeling (to account for few numbers of replicates and large dynamic expression range). 

# Distributions for count data

Explain when to use common statistical models for counts:

Poisson - normal distribution, large sample number (microarray data that has dynamic range limited maximum due to when the probes max out, and therefore uses the Poisson). Why can we use Poisson for microarray when also there are small sample sizes?

NB - non-normal distribution, low number of biological reps. No max for dynamic range. Don't use technical replicates - biological replicates much more useful (cite paper? or later in md)


# DESeq2 workflow

The DESeq2 workflow is shown below in green. After obtaining the counts associated with each gene, DESeq2 normalizes the count values to account for differences in library sizes and RNA composition between samples. Then, QC is performed at the gene and sample level prior to performing the differential expression analysis.

<img src="../img/deseq_workflow_full.png" width="200">

## Normalization

The first step in the workflow is count normalization, which is necessary to make accurate comparisons of gene expression between samples. The raw counts, or number of reads aligning to each gene, need to be normalized to account for differences in library depth between samples when performing differential expression analyses.

<img src="../img/deseq_workflow_normalization.png" width="200">

While normalization is necessary for differential expression analyses, it is also necessary whenever exploring or comparing counts between or within samples. 

Different types of normalization methods exist, and a few of the most common methods include:
 
 - **normalization for library size:** necessary for comparison of the same gene between samples
 
 	<img src="../img/sequencing_depth.png" width="400">
 
 - **normalization for gene length:** necessary for comparison of different genes within a sample
 
 	<img src="../img/length_of_gene.png" width="400">
 
 - **normalization for RNA composition:** recommended for comparison between samples (particularly important when performing differential expression analyses)
 
 	>"A few highly and differentially expressed genes may have strong influence on the total read count, causing the ratio of total read counts not to be a good estimate for the ratio of expected counts (for all genes)"[[1](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)]
 
 
### Common normalization measures

Several common normalization measures exist to account for these differences:

- **CPM (counts per million):** counts scaled by total number of reads. This measure accounts for sequencing depth only.
- **TPM (transcripts per kilobase million):** counts per length of transcript (kb) per million reads mapped. This measure accounts for both sequencing depth and gene length.
- **RPKM/FPKM (reads/fragments per kilobase of exon per million reads/fragments mapped):** similar to TPM, as this measure accounts for both sequencing depth and gene length as well; however, it is **not recommended**.
- **Tool-specific metrics for normalization:** 
	- DESeq2 uses a median of ratios method, which accounts for sequencing depth and RNA composition [[1](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)]. 
	- EdgeR uses a trimmed mean of M values (TMM) method that accounts for sequencing depth, RNA composition, and gene length [[2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)]

### RPKM/FPKM (not recommended)
While TPM and RPKM/FPKM normalization methods both account for sequencing depth and gene length, RPKM/FPKM measures are not recommended. **The reason  is that the normalized count values output by the RPKM/FPKM method are not comparable between samples.** 

Using RPKM/FPKM normalization, the total number of RPKM/FPKM normalized counts for each sample will be different. Therefore, you cannot compare the normalized counts for each gene equally between samples. 

**RPKM-normalized counts table**

| gene | sampleA | sampleB |
| ----- |:-----:|:-----:|
| MOV10 | 5.5 | 5.5 |
| ABCD | 73.4 | 21.8 |
| ... | ... | ... |
|Total RPKM-normalized counts | 1,000,000 | 1,500,000 |

SampleA has a greater proportion of counts associated with MOV10 (5.5/1,000,000) than does sampleB (5.5/1,500,000) even though the RPKM count values are the same. Therefore, we cannot directly compare the counts for MOV10 (or any other gene) between sample1 and sample2 because the total number of normalized counts are different between samples. 

### TPM (recommended)
In contrast to RPKM/FPKM, TPM-normalized counts normalize for both sequencing depth and gene length, but have the same total TPM-normalized counts per sample. Therefore, the normalized count values are comparable both between and within samples.

> *NOTE:* [This video by StatQuest](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) shows in more detail why TPM should be used in place of RPKM/FPKM if needing to normalize for sequencing depth and gene length.

### DESeq2-normalized counts - Median of ratios method
Since tools for differential expression analysis are comparing the counts between sample groups for the same gene, gene length does not need to be accounted for by the tool. However, **sequencing depth** and **RNA composition** do need to be taken into account.

To normalize for sequencing depth and RNA composition, DESeq2 uses the median of ratios method, which performs the following steps when you run the tool:

**Step 1: creates a pseudo-reference sample (row-wise geometric mean)**

For each gene, a pseudo-reference sample is created that is equal to the geometric mean across all samples.

| gene | sampleA | sampleB | pseudo-reference sample  |
| ----- |:-----:|:-----:|:-----:|
| EF2A | 1489 | 906 | sqrt(1489 * 906) = **1161.5** |
| ABCD | 22 | 13 | sqrt(24 * 13) = **17.7** |
| ... | ... | ... | ... |

**Step 2: calculates ratio of each sample to the reference**

For every gene in a sample, the ratios (sample1/ref) are calculated (as shown below). This is performed for each sample in the dataset. Since the majority of genes are not differentially expressed, the majority of genes in each sample should have similar ratios within the sample.

| gene | sampleA | sampleB | pseudo-reference sample  | ratio sample1/ref | ratio sample2/ref |
| ----- |:-----:|:-----:|:-----:| :-----: | :-----: |
| EF2A | 1489 | 906 | 1161.5 | 1489/1161.5 = **1.28** | 906/1161.5 = **0.78** |
| ABCD | 22 | 13 | 16.9 | 22/16.9 = **1.30** | 13/16.9 = **0.77** |
| MEF3 | 793 | 410 | 570.2 | 793/570.2 = **1.39** | 410/570.2 = **0.72**
| BBC1 | 76 | 42 | 56.5 | 76/56.5 = **1.35** | 42/56.5 = **0.74**
| MOV10 | 521 | 1196 | 883.7 | 521/883.7 = **0.590** | 1196/883.7 = **1.35** |
| ... | ... | ... | ... |

**Step 3: takes sample's median value as that sample's normalization factor**

The median value of all ratios is taken as the normalization factor (size factor) for that sample, as calculated for SampleA below. Notice that the differentially expressed genes should not affect the median value:

`normalization_factor_sample1 <- median(c(0.59, 1.28, 1.3, 1.35, 1.39))`

`normalization_factor_sample2 <- median(c(0.72, 0.74, 0.77, 0.78, 1.35))`
 
The figure below illustrates the median value for the distribution of all gene ratios for a single sample (frequency is on the y-axis).

<img src="../img/deseq_median_of_ratios.png" width="400">

The median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample (large outlier genes will not represent the median ratio values). **This method is robust to imbalance in up-/down-regulation and large numbers of differentially expressed genes.**

**Step 4: divide each raw count value in sample by that sample's normalization factor to generate normalized count values**

For example, if median ratio for SampleA was 1.3 and the median ratio for SampleB was 0.77, you could calculate normalized counts as follows:

SampleA median ratio = 1.3

SampleB median ratio = 0.77

**Raw Counts**

| gene | sampleA | sampleB |  
| ----- |:-----:|:-----:|
| EF2A | 1489 | 906 | 
| ABCD | 22 | 13 | 
| ... | ... | ... | 

**Normalized Counts**

| gene | sampleA | sampleB |
| ----- |:-----:|:-----:|
| EF2A | 1489 / 1.3 = **1145.39** | 906 / 0.77 = **1176.62** | 
| ABCD | 22 / 1.3 = **16.92** | 13 / 0.77 = **16.88** | 
| ... | ... | ... | 

***
**Exercise**

Determine the normalized counts for your gene of interest, PD1, given the raw counts and size factors below. 

NOTE: You will need to run the code below to generate the raw counts dataframe (PD1) and the size factor vector (size_factors), then use these objects to determine the normalized counts values:

```r

# Raw counts for PD1
PD1 <- c(21, 58, 17, 97, 83, 10)
names(PD1) <- paste("Sample", 1:6)
PD1 <- data.frame(PD1)
PD1 <- t(PD1)

# Size factors for each sample
size_factors <- c(1.32, 0.70, 1.04, 1.27, 1.11, 0.85)

```

***

## Quality Control

The next step in the DESeq2 workflow is QC, which includes sample-level and gene-level steps.

<img src="../img/deseq_workflow_qc.png" width="200">

### Sample-level QC

A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design? Log2-transformed normalized counts are used to assess similarity between samples using Principal Component Analysis (PCA) and hierarchical clustering.

<img src="../img/sample_qc.png" width="600">


#### PCA

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). Details regarding PCA are given below (based on [materials from StatQuest](https://www.youtube.com/watch?v=_UVHneBUBW0), and if you would like a more thorough description, we encourage you to explore [StatQuest's video](https://www.youtube.com/watch?v=_UVHneBUBW0). 

If you had two samples and wanted to plot the counts of one sample versus another, you could plot the counts of one sample on the x-axis and the other sample on the y-axis as shown below:

<img src="../img/PCA_2sample_genes.png" width="600">

You could draw a line through the data in the direction representing the most variation, which is on the diagonal in this example. The maximum variation in the data is between the two endpoints of this line.  

We also see the genes vary somewhat above and below the line. We could draw another line through the data representing the second most amount of variation in the data. 


<img src="../img/PCA_2sample_variation1.png" width="600">

The genes near the ends of the line, which would include those genes with the highest variation between samples (high expression in one sample and low expression in the other), have the greatest influence on the direction of the line. 

<img src="../img/PCA_2sample_variation2.png" width="600">

For example, a small change in the value of *Gene C* would greatly change the direction of the line, whereas a small change in *Gene A* or *Gene D* would have little affect.

<img src="../img/PCA_2sample_variation3.png" width="900">

We could just rotate the entire plot and view the lines representing the variation as left-to-right and up-and-down. We see most of the variation in the data is left-to-right; this is and the second most variation in the data is up-and-down. These axes that represent the variation are "Principal Components", with PC1 representing the most variation in the data and PC2 representing the second most variation in the data. 

If we had three samples, then we would have an extra direction in which we could have variation. Therefore, if we have *N* samples we would have *N*-directions of variation or principal components.

<img src="../img/PCA_2sample_rotate.png" width="300">

We could give quantitative scores to genes based on how much they influence PC1 and PC2. Genes with little influence would get scores near zero, while genes with more influence would receive larger scores. Genes on opposite ends of the lines have a large influence, so they would receive large scores, but with opposite signs.

<img src="../img/PCA_2sample_influence.png" width="600">

To generate a score per sample, we combine the read counts for all genes. To calculate the scores, we do the following:
	
	Sample1 PC1 score = (read count * influence) + ... for all genes
	
Using the counts in the table for each gene (assuming we had only 4 genes total) we could calculate PC1 and PC2 values for each sample as follows:

	Sample1 PC1 score = (4 * -2) + (1 * -10) + (8 * 8) + (5 * 1) = 51
	Sample1 PC2 score = (4 * 0.5) + (1 * 1) + (8 * -5) + (5 * 6) = -7
	
	Sample2 PC1 score = (5 * -2) + (4 * -10) + (8 * 8) + (7 * 1) = 21
	Sample2 PC2 score = (5 * 0.5) + (4 * 1) + (8 * -5) + (7 * 6) = 8.5
	
The scores would then be plotted to examine whether the samples exhibit similar variation across all genes:

<img src="../img/PCA_samples.png" width="600">

Since genes with the greatest variation between samples will have the greatest influence on the principal components, we hope our condition of interest explains this variation (e.g. high counts in one condition and low counts in the other). With PC1 representing the most variation in the data and PC2 representing the second most variation in the data, we can visualize how similar the variation of genes is between samples.

#### Hierarchical Clustering Heatmap

Similar to PCA, hierarchical clustering is another, complementary method for identifying strong patterns in a dataset and potential outliers. The heatmap displays the correlation or distances for all pairwise combinations of samples. Since the majority of genes are not differentially expressed, samples generally have high correlations with each other (values higher than 0.80). Samples below 0.80 may indicate an outlier in your data and/or sample contamination.

The hierarchical tree can indicate which samples are more similar to each other based on the normalized gene expression values. The color blocks indicate substructure in the data, and you would expect to see your replicates separate together as a block for each sample group.

<img src="../img/heatmap_example.png" width="500">

### Gene-level QC

Prior to differential expression analysis it is beneficial to omit genes that have little or no chance of being detected as differentially expressed. This will increase the power to detect differentially expressed genes. The genes omitted fall into three categories:

- Genes with zero counts in all samples
- Genes with an extreme count outlier
- Genes with a low mean normalized count across all samples

<img src="../img/gene_filtering.png" width="600">

**DESq2 will perform this filtering by default; however other DE tools, such as EdgeR will not.** It is important to understand what filtering is performed by your tool of choice to know if you need to perform any additional filtering prior to the differential expression analysis.

## Differential expression analysis with DESeq2

The final step in the differential expression analysis is the actual fitting of the raw counts to the statistical model and testing for differentially expressed genes. Essentially we want to determine whether the mean expression levels of two different samplegroups are significantly different.

<img src="../img/de_theory.png" width="600">


To determine differentially expressed genes, we are going to use the DESeq2 tool.   This tool builds on good ideas for dispersion estimation and use of Generalized Linear Models from the DSS and edgeR methods. The [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) was published in 2014, but the package is continually updated and available for use in R through Bioconductor.

Differential expression analysis with DESeq2 requires multiple steps, as displayed below.

<img src="../img/DESeq2_workflow.png" width="500">

### Estimate size factors

As you have probably noticed, this is exactly what we did to normalize our counts. DESeq2 will automatically estimate the size factors when performing the differential expression analysis if you haven't already done so. If you have already generated the size factors, then DESeq2 will use these values.

### Estimate variation / dispersion

To accurately model our sequencing counts, we need to generate accurate estimates of within-group variation (variation between replicates of the same samplegroup) for each gene. With only a few (3-6) replicates per group, the estimates of variation for each gene are often unreliable. Therefore, DESeq2 shares information across genes to generate more accurate estimates of variation based on the expression level of the gene using a method called 'shrinkage'. DESeq2 assumes that genes with similar expression levels have similar dispersion or variation of expression. DESeq2 generates more accurate measures of dispersion using the following steps:

**Step 1: Estimate the dispersion for each gene separately**

To model the dispersion based on expression level (mean normalized counts of replicates), the dispersion for each gene is estimated using maximum likelihood estimation. Given the normalized count values of the replicates, the most likely estimate of dispersion is calculated.

**Step 2: Fit a curve to the the gene estimates given expression strength**

The idea behind fitting a curve to the data is that different genes will have different scales of biological variability, but, over all genes, there will be a distribution of reasonable estimates of dispersion. 

**Step 3: Shrink gene-wise dispersion estimates toward the values predicted by the curve**

The curve allows for more accurate identification of differentially expressed genes when sample sizes are small, and the strength of the shrinkage for each gene depends on :
	
- how close gene dispersions are from the curve
- sample size (more samples = less shrinkage)


This shrinkage method is particularly important to reduce false positives in the differential expression analysis. Genes with extremely low levels of variation are shrunken towards the curve, and the more accurate, higher dispersion values are output for differential expression testing. 

This curve is displayed as a red line in the figure below, which plots the estimate for the expected dispersion value for genes of a given expression strength. Each black dot is a gene with an associated mean expression level and dispersion estimate.

Genes with extremely high dispersion values are not shrunken toward the curve due to the likelihood that the gene does not follow the modeling assumptions and has higher variability than others for biological or technical reasons. Shrinking the values toward the curve could result in false positives, so these values are not shrunken. These genes are shown surrounded by blue circles below. 

<img src="../img/deseq_dispersion.png" width="500">

> **NOTE:** This is a good plot to ensure your data is a good fit for the DESeq2 model. You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing expression levels. If you see a cloud or different shapes, then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.

### Generalized Linear Model fit for each gene

As discussed earlier, the count data generated by RNA-Seq exhibits overdispersion and the statistical distribution used to model the counts needs to account for this over dispersion. DESeq2 uses a negative binomial distribution to model the RNA-Seq counts using the equation below:

 <img src="../img/NB_model_formula.png" width="500">
 
DESq2 will use this formula to create the model for each gene, but what we really want to know is the log2 foldchanges between conditions. The log2 foldchanges can be estimated by replacing the normalized counts in the model using the formula:

 <img src="../img/NB_model_formula_betas.png" width="600">

By fitting the model, DESeq2 will determine the estimates for the log2 foldchanges between conditions. 

For example, if we had our samples divided into three conditions, control (ctrl), overexpression (oe) and knockdown (kd), our model would be:

 <img src="../img/NB_model_formula_betas_example.png" width="570">

In DESeq2 and most other DE tools, you will assign your samples to specific conditions using a 'model matrix' or 'design matrix'. For example, we can assign our samples to conditions in the model matrix using binary (0,1) notation:

|  | ctrl | oe | kd |
| ----- |:-----:|:-----:|:-----:|
| sample1 | 1 | 0 | 0 |
| sample2 | 1 | 0 | 0 |
| sample3 | 0 | 1 | 0 |
| sample4 | 0 | 1 | 0 |
| sample5 | 0 | 0 | 1 |
| sample6 | 0 | 0 | 1 |

Samples 1 and 2 are controls, samples 3 and 4 are overexpression, and samples 5 and 6 are knockdown. This information is utilized to inform the model about which replicates should be used to estimate each of the log2 foldchanges. For example, sample1 and sample2 should be used to estimate the log2 foldchanges for the control group (relative to mean expression of all groups), and the model formula for these samples would be:

<img src="../img/NB_model_formula_betas_example3.png" width="400">

We will discuss later how to obtain log2 foldchanges relative to other sample groups instead of to the mean expression of all groups.

##Fold change shrinkage (Fisher info: degrees of freedom (number of samples, number of betas), estimated mean counts, dispersion estimate); beta prior (briefly - if very few reps and high variation within gene); Wald testing; LRT testing; exercises
 

