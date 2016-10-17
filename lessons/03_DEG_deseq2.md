---
title: "Gene-level differential expression analysis using DESeq2"
author: "Meeta Mistry"
date: "October 17, 2016"
---

Approximate time: 45 minutes

## Learning Objectives 

* Understanding the different components of differential expression analysis in the context of DESeq2
* Exploring different objects in DESeq2 
* Summarizing and filtering results to find significant DEGs


## DESeq2: Differential expression analysis

### Getting setup

Let's get started by opening RStudio and opening up the project that we created last lesson. 

1. Go to the File menu and select 'Open project ...'
2. Navigate to `~/Desktop/DEanalysis/` and double click on the `DEanalysis.Rproj` file

You should see your environment become populated with all of the variables created last lesson. The only thing that we will need to do is reload the required libraries:

```
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```


### Running DESeq2

To run the differential expression pipeline on the raw counts in DESeq2, we use a **single call to the function `DESeq()`**. The required input is the `DESeqDataSet` object that we created in the last lesson. By re-assigning the results of the function back to the same variable name, we can continue to fill in the `slots` of our `DESeqDataSet` object.

	##Run analysis
	dds <- DESeq(dds)
 
This function will print out a message for the various steps it performs: 

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
``` 

<img src="../img/slide16+33_DGE.png" width="600">


**Everything from normalization to linear modeling was carried out by the use of a single function!** The results of each step were inserted into the object that you initialized.

![deseq1](../img/deseq_obj2.png)


> *NOTE:* There are individual functions available in DESeq2 that would allow us to carry out each step in the workflow in a step-wise manner, rather than a single call. We demonstrated one example when generating size factors to create a normalized matrix. By calling `DESeq()`, the individual functions for each step are run for you.


## Normalization

To normalize the count data DESeq2 calculates size factors for each sample, using the *median of ratios method*. Let's take a quick look at size factor values we have for each sample:

```
> sizeFactors(dds)
Mov10_kd_2 Mov10_kd_3 Mov10_oe_1 Mov10_oe_2 Mov10_oe_3 Irrel_kd_1 Irrel_kd_2 Irrel_kd_3 
 1.5646728  0.9351760  1.2016082  1.1205912  0.6534987  1.1224020  0.9625632  0.7477715 
 
```
 
These numbers should be identical to those we generated initially when we had run the function `estimateSizeFactors(dds)`. Take a look at the total number of reads for each sample using `colSums(counts(dds))`. *How do the numbers correlate with the size factor?*

> *NOTE:* it can be advantageous to calculate gene-specific normalization factors (size factors) to account for further sources of technical biases such as differing dependence on GC content, gene length or the like, and these can be supplied to DESeq2 instead of using the median of ratios method.


## Dispersion estimates

In our model, the **within group variability** is accounted for using the dispersion parameter. Dispersion estimates are computed **per gene**, because different genes naturally have a different scale of biological variability. DESeq2 does a first pass estimate on dispersion for each gene (using maximum-likelihood estimate), but with such small sample sizes we will make very bad estimates of gene-wise dispersion unless we **share information across genes**. The next step is therefore taking information from all gene dispersion estimates to shrink them to more reasonable values.

Let's take a look at the dispersion estimates for our data:

	# Plot dispersion estimates
	plotDispEsts(dds)
	

<img src="../img/plotDispersion.png"">
 

The black dots are the original estimates for each gene. The red smooth curve provides an accurate estimate for the expected dispersion value for genes of a given expression strength. The blue dots represent shrunken estimates. The circles indicate outliers, where we don't perform shrinkage. 

We use an empirical Bayes approach which lets the strength of shrinkage depend (i) on an estimate of how close true dispersion values tend to be to the fit and (ii) on the degrees of freedom. **Since we have a small sample size, for many genes we see quite a bit of shrinkage.**


## Identifying gene expression changes

We have three sample classes so we can make three possible pairwise comparisons:

1. Control vs. Mov10 overexpression
2. Control vs. Mov10 knockdown
3. Mov10 knockdown vs. Mov10 overexpression

**We are really only interested in #1 and #2 from above**. Using the design formula we provided `~sampletype`, DESeq 2 internally created the following design matrix:

```
   	      Intercept	sampletypecontrol sampletypeMOV10_knockdown	sampletypeMOV10_overexpression
Mov10_kd_2	 1		0		1		0
Mov10_kd_3	 1		0		1		0
Mov10_oe_1   1		0		0		1
Mov10_oe_2   1		0		0		1
Mov10_oe_3   1		0		0		1
Irrel_kd_1	 1		1		0		0
Irrel_kd_2	 1		1		0		0
Irrel_kd_3	 1		1		0		0	

```
This design matrix is now used to setup the contrasts to request the comparisons we want to make.


### Hypothesis testing: Wald test

To build a results table, we use the `results()` function on the `dds` object. Additionally we need to specify **which comparisons we are interested in** looking at. 

The comparisons are provided to DESeq2 in the form of **contrasts**, in one of three different ways. In this lesson we will demonstrate the method that is most intuitive. By providing contrasts we are telling DESeq2 **which coefficients to use for the hypothesis testing** procedure; this also corresponds to the headers in your design matrix. To find out how the coefficients are named we can use the `resultsNames()` function:

	# Find names of coefficients
	resultsNames(dds)

To specify the specific contrasts we are interested in, we need to provide the column names from the coefficients table as a list of 2 character vectors:

	## Define contrasts
	contrast_oe <- list( "sampletypeMOV10_overexpression", "sampletypecontrol")

**The order of the names, determines the direction of fold change that is reported.** The name provided in the second element is the level that is used to baseline. So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control. Pass the contrast vector as an argument to the `results()` function:

	# Extract results table
	res_tableOE <- results(dds, contrast=contrast_oe)


Let's take a look at what information is stored in the results:

	head(res_tableOE)

```
log2 fold change (MAP): sampletype MOV10_overexpression vs control 
Wald test p-value: sampletype MOV10_overexpression vs control 
DataFrame with 6 rows and 6 columns
               baseMean log2FoldChange      lfcSE       stat    pvalue       padj
              <numeric>      <numeric>  <numeric>  <numeric> <numeric>  <numeric>
1/2-SBSRNA4  45.6520399     0.26976764 0.18775752  1.4367874 0.1507784 0.25242910
A1BG         61.0931017     0.20999700 0.17315013  1.2128030 0.2252051 0.34444163
A1BG-AS1    175.6658069    -0.05197768 0.12366259 -0.4203185 0.6742528 0.77216278
A1CF          0.2376919     0.02237286 0.04577046  0.4888056 0.6249793         NA
A2LD1        89.6179845     0.34598540 0.15901426  2.1758136 0.0295692 0.06725157
A2M           5.8600841    -0.27850841 0.18051805 -1.5428286 0.1228724 0.21489067
```
> *NOTE:* The results table looks very much like a data frame and in many ways it can be treated like one (i.e when accessing/subsetting data). However, it is important to recognize that it is actually stored in a `DESeqResults` object. When we start visualizing our data, this information will be helpful. 


Let's go through some of the columns in the results table to get a better idea of what we are looking at. To extract information regarding the meaning of each column we can use `mcols()`:

	mcols(res_tableOE, use.names=T)

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 

***

**Exercise**

1. Create a contrasts vector called `contrast_kd` for the Mov10_knockdown comparison to control.
2. Use that contrasts vector to extract a results table and store that to a variable called `res_tableKD`.  
3. Create a contrasts vector for the Mov10_overexpression comparison to *all other samples*.

*** 


### Summarizing results and identifying DEGs

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results at a given FDR threshold. 

	## Summarize results
	summary(res_tableOE)
	

```  
out of 19748 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 3657, 19% 
LFC < 0 (down)   : 3897, 20% 
outliers [1]     : 0, 0% 
low counts [2]   : 3912, 20% 
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

In addition to the number of genes up- and down-regulated at and FDR < 0.1, the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count (which in our case is < 4).

The default FDR threshold is set to `alpha = 0.1`, which is quite liberal. Let's try changing that to `0.05` -- *how many genes are we left with*?

The FDR threshold on it's own doesn't appear to be reducing the number of significant genes. With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also add a fold change threshold. The `summary()` function doesn't have an argument for fold change threshold, but instead we can use the base R function `subset()`.

Let's first create variables that contain our threshold criteria:

	### Set thresholds
	padj.cutoff <- 0.05
	lfc.cutoff <- 1

The `lfc.cutoff` is set to 1; remember that we are working with log2 fold changes so this translates to an actual fold change of 2 which is pretty reasonable. Now let's setup our **`subset()` function**. Start building from the inside out:

	subset(res_tableOE)

We need to add our selection criteria. The first is our FDR threshold:

	subset(res_tableOE, padj < padj.cutoff)

Now let's add in the log2 fold change criteria. Because we want both up- and down-regulated genes we will use the absolute value of the fold change using the `abs(log2FoldChange)` function:

	subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

Now, finally we will put all of that inside the `summary()` function. This is a fast way of getting overall statistics and deciding whether our threshold is still too liberal or perhaps overly stringent.

	summary(subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff), alpha =0.05)


Does this reduce our results? How many genes are up-regulated and down-regulated at this new threshold?

We should have a total of 884 genes (682 up-regulated and 202 down-regulated) that are significantly differentially expressed. To denote these genes as significant we can add a column in our results table. The column will be a logical vector, where `TRUE` means the gene passes our threshold and `FALSE` means it fails.

	# Add a threshold vector
	threshold <- res_tableOE$padj < padj.cutoff & 
                   abs(res_tableOE$log2FoldChange) > lfc.cutoff
                   
To add this vector to our results table we can use the `$` notation to create the column on the left hand side of the assignment operator, and the assign the vector to it:

	res_tableOE$threshold <- threshold                

Now we can easily check how many genes are significant by using the `which()` function:

	length(which(res_tableOE$threshold))

***

**Exercise**

1. Explore the results table summary for the **Mov10_knockdown comparison to control**. How many genes are differentially expressed using the default thresholds?
2. Using the same thresholds as above (`padj.cutoff < 0.05` and `lfc.cutoff = 1`), report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.
3. Add a new column called `threshold` to the `res_tableKD` which contains a logical vector denoting genes as being differentially expressed or not.

*** 

> **NOTE: on p-values set to NA**
> > 
> 1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
> 2. If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance. 
> 3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA. 




---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*
