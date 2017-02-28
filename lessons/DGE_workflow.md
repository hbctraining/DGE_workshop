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

The first step in the workflow is normalization, which is necessary to make accurate comparisons of gene expression between samples. The raw counts, or number of reads aligning to each gene, need to be normalized to account for differences in library depth between samples, such as when performing differential expression analyses. If making comparisons between genes within samples, gene length should also be taken into account. 

Different types of normalization methods exist, and a few of the most common methods include:
 
 - **normalization for library size:** comparison between samples
 
 	<img src="../img/sequencing_depth.png" width="400">
 
 - **normalization for gene length:** comparison within samples
 
 	<img src="../img/length_of_gene.png" width="400">
 
 - **normalization for RNA composition:** comparison between samples (particularly important for differential expression analyses)
 
 	>"A few highly and differentially expressed genes may have strong influence on the total read count, causing the ratio of total read counts not to be a good estimate for the ratio of expected counts (for all genes)"[[1](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)]
 
 
### Common normalization measures

Several common normalization measures exist to account for these differences:

- **CPM (counts per million):** counts scaled by total number of reads. This measure accounts for sequencing depth only.
- **TPM (transcripts per kilobase million):** counts per length of transcript (kb) per million reads mapped. This measure accounts for both sequencing depth and gene length.
- **RPKM/FPKM (reads/fragments per kilobase of exon per million reads/fragments mapped):** similar to TPM. This measure accounts for both sequencing depth and gene length as well; however, it is **not recommended**.
- **Tool-specific metrics for normalization:** 
	- DESeq2 uses a median of ratios method, which accounts for sequencing depth and RNA composition[[1]()]. 
	- EdgeR uses a trimmed mean of M values (TMM) method that accounts for sequencing depth, RNA composition, and gene length [[2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)]

#### RPKM/FPKM (not recommended)
While TPM and RPKM/FPKM normalization methods both account for sequencing depth and gene length, RPKM/FPKM measures are not recommended. **The reason  is that the normalized count values output by the RPKM/FPKM method are not comparable between samples.** 

If you sum the total number of RPKM/FPKM normalized counts for each sample, the total number will be different between samples. Therefore, you cannot compare the normalized counts for each gene equally between samples. 


**RPKM-normalized counts table**

| gene | sample1 | sample2 |
| ----- |:-----:|:-----:|
| MOV10 | 5.5 | 5.5 |
| ABCD | 73.4 | 21.8 |
| ... | ... | ... |
|Total RPKM-normalized counts | 1,000,000 | 1,500,000 |

For example, using the table above, the MOV10 gene has 5.5 RPKM-normalized counts in sample1 out of 1,000,000 total RPKM-normalized counts in the sample, while sample2 also has 5.5 RPKM-normalized counts associated with MOV10 out of 1,500,000 total RPKM-normalized counts for sample2. Therefore, we cannot directly compare the counts for MOV10 (or any other gene) between sample1 and sample2 because the total number of normalized counts are different between samples. 

Sample1 has a greater proportion of counts associated with MOV10 (5.5/1,000,000) than does sample2 (5.5/1,500,000) even though the RPKM count values are the same.

> *NOTE:* [This video by StatQuest](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) shows in more detail why TPM should be used in place of RPKM/FPKM if needing to normalize for sequencing depth and gene length.

#### TPM (recommended)
In contrast to RPKM/FPKM, TPM-normalized counts normalize for both sequencing depth and gene length, but have the same total TPM-normalized counts per sample. Therefore, the normalized count values are comparable both between and within samples.

#### DESeq2-normalized counts - Median of ratios method
Since tools for differential expression analysis are comparing the counts between sample groups for the same gene, gene length does not need to be accounted for by the tool. However, sequencing depth and RNA composition do need to be included in the normalization.

To normalize for sequencing depth and RNA composition, DESeq2 uses the median of ratios method, which performs the following steps when you run the tool:

**Step 1: creates a pseudo-reference sample (row-wise geometric mean)**

For each gene, a pseudo-reference sample is created that is equal to the geometric mean across all samples.

| gene | sample1 | sample2 | pseudo-reference sample  |
| ----- |:-----:|:-----:|:-----:|
| MOV10 | 1489 | 906 | sqrt(1489 * 906) = **1161.5** |
| ABCD | 24 | 13 | sqrt(24 * 13) = **17.7** |
| ... | ... | ... | ... |

**Step 2: calculates ratio of each sample to the reference**

For every gene in a sample, the ratios (sample1/ref) are calculated (as shown below). This is performed for each sample in the dataset.

| gene | sample1 | sample2 | pseudo-reference sample  | ratio sample1/ref | ratio sample2/ref |
| ----- |:-----:|:-----:|:-----:| :-----: | :-----: |
| EF2A | 1489 | 906 | 1161.5 | 1489/1161.5 = **1.28** | 906/1161.5 = **0.78** |
| BBC1 | 22 | 13 | 16.9 | 22/16.9 = **1.30** | 13/16.9 = **0.77** |
| ... | ... | ... | ... |

**Step 3: takes sample's median value as that sample's normalization factor**

 The median value of all ratios is taken as the normalization factor (size factor) for that sample, as depicted in the figure below. The figure shows a histogram with the size factors versus frequency of genes with those ratios. 

<img src="../img/deseq_median_of_ratios.png" width="400">

The median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample (large outlier genes will not represent the median ratio values). This method is robust to imbalance in up-/down-regulation and large numbers of differentially expressed genes.

## QC (sample-level and gene-level)

A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). The explanation of PCA below utilizes materials available from StatQuest, and if you would like a more thorough description, we encourage you to explore the video. If you had two samples and wanted plot the counts of one sample versus another, you could do the following:

<img src="../img/PCA_2sample_genes.png" width="600">

You could draw a line through the data representing the most variation in the data, which is on the diagonal in this dataset. The maximum variation in the data is between the two endpoints of this line.  We also see the genes vary somewhat above or below the line. We could draw another line through the data representing the second most amount of variation in the data. 


<img src="../img/PCA_2sample_variation1.png" width="600">

The genes near the ends of the line, which would inlude those genes with the highest variation between samples (high expression in one sample and low expression in the other), have the greatest influence on the direction of the line. 

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

## Variation / Dispersion

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

> **NOTE:** This is a good plot to ensure your data is a good fit for the DESeq2 model. You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing expression levels. If you see a cloud of data or different shapes, then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.
