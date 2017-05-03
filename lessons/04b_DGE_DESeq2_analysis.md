







### Multiple test correction

If we used the p-value directly from the Wald test with a significance cut-off of 0.05 (α = 0.05), then 5% of all genes would be called as differentially expressed (i.e. 5% False positive genes). The more genes we test, the more 'false positives' we discover. DESeq2 helps reduce the number of genes tested by removing those genes unlikely to be significantly DE, such as those with low number of counts and outlier samples. However, we still need to correct for multiple testing, and there are a few common approaches:

- **Bonferroni:** Reject any hypothesis with p-value ≤ α/m. **This is a very conservative approach with a high probability of false negatives.**
- **FDR / Benjamini-Hochberg:** Rank j / m multiplied by the FDR levels. This approach is designed to control the proportion of false positives among the set of rejected hypotheses
- **Q-value:** The minimum FDR that can be attained when calling that feature significant. For example, if gene X has a q-value of 0.013 it means that 1.3% of genes that show p-values at least as small as gene X are false positives

In DESeq2, the p-values attained by the Wald test are corrected for multiple testing using the Benjamin and Hochberg method. The p-adjusted values should be used to determine significant genes. The significant genes can be output for visualization and/or functional analysis.

## Differential expression analysis of Mov10 dataset 

Let's put the theory into practice by performing differential gene expression analysis on the Mov10 dataset using DESeq2. We will be using several different R packages to perform the analysis, so we will start by loading these libraries:

### Loading libraries

```r
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```



>**NOTE:** At this step, if you expected a number of genes to have log2 fold changes far outside the normal range for your dataset, you could turn off the beta prior (as discussed previously), which would turn off the LFC shrinkage. You can turn off the beta prior using `DESeq(dds, betaPrior=FALSE)`.



## Normalization

To normalize the count data DESeq2 calculates size factors for each sample, using the *median of ratios method* discussed previously. Let's take a quick look at size factor values we have for each sample:

```
sizeFactors(dds)

Mov10_kd_2 Mov10_kd_3 Mov10_oe_1 Mov10_oe_2 Mov10_oe_3 Irrel_kd_1 Irrel_kd_2 Irrel_kd_3 
 1.5646728  0.9351760  1.2016082  1.1205912  0.6534987  1.1224020  0.9625632  0.7477715  
```
 
These numbers should be identical to those we generated initially when we had run the function `estimateSizeFactors(dds)`. Take a look at the total number of reads for each sample using:

```r
colSums(counts(dds))
```

*How do the numbers correlate with the size factor?*

Now take a look at the total depth after normalization using:

```r
colSums(counts(dds, normalized=T))
```
How do the values across samples compare with the total counts taken for each sample?

> **NOTE:** it can be advantageous to calculate gene-specific normalization factors (size factors) to account for further sources of technical biases such as differing dependence on GC content, gene length or the like, and these can be supplied to DESeq2 instead of using the median of ratios method.

## Dispersion estimates

In our model, the **within group variability** is accounted for using the dispersion parameter. Dispersion estimates are computed **per gene**, because different genes naturally have a different scale of biological variability. The next step is therefore taking information from all gene dispersion estimates to shrink them to more reasonable values.

Let's take a look at the dispersion estimates for our data:

```r
# Plot dispersion estimates
plotDispEsts(dds)
```

<img src="../img/plotDispersion.png">
 

The black dots are the original estimates for each gene. The red smooth curve provides an accurate estimate for the expected dispersion value for genes of a given expression strength. The blue dots represent shrunken estimates. The circles indicate outliers, where we don't perform shrinkage. Remember that the strength of shrinkage depend on:
	
- how close gene dispersions are from the curve
- sample size (more samples = less shrinkage)
	
**Since we have a small sample size, for many genes we see quite a bit of shrinkage.**

This is often a plot we check to make sure that our data are a good fit for the DESeq2 model based on how well the dispersion estimates follow the curve. **Do you think our data are a good fit?**


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
This design matrix is now used to setup the contrasts to request the comparisons we want to make. This information is utilized to inform the model about which replicates should be used to estimate the **log2 foldchanges (LFC)**.



### Hypothesis testing: Wald test

DESeq2 performs a hypothesis test for all possible pairwise comparisons. In order for us to retrieve the results for a specific pair of sample classes we need to provide this information to DESeq2 in the form of contrasts. While contrasts can be provided a few different ways, we will use the `list()` method:

We need to use the coefficient names to specify our comparisons, these correspond to the headers in your design matrix. To find out how the coefficients are named we can use the resultsNames() function:


```r
# Find names of coefficients
resultsNames(dds)
```

To specify the specific contrasts, we need to provide the column names from the coefficients table as a list of 2 character vectors:

```r
## Define contrasts
contrast_oe <- list( "sampletypeMOV10_overexpression", "sampletypecontrol")
```

**The order of the names, determines the direction of fold change that is reported.** The name provided in the second element is the level that is used to baseline. So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control. Pass the contrast vector as an argument to the `results()` function:

```r
# Extract results table
res_tableOE <- results(dds, contrast=contrast_oe)
```

>**NOTE:** We could have specified the contrast in the `results()` argument:
`results(dds, contrast = c("sampletype", "MOV10_overexpression", "control	normal"))`.

This will build a results table containing Wald test statistics for the comparison we are interested in. Let's take a look at what information is stored in the results:

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

The results table looks very much like a dataframe and in many ways it can be treated like one (i.e when accessing/subsetting data). However, it is important to recognize that it is actually stored in a `DESeqResults` object. When we start visualizing our data, this information will be helpful. 

```r
class(res_tableOE)
```

Let's go through some of the columns in the results table to get a better idea of what we are looking at. To extract information regarding the meaning of each column we can use `mcols()`:

```r
mcols(res_tableOE, use.names=T)
```

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 
Now that we have results for the overexpression results, let's do the same for the **Control vs. Knockdown samples**. The first thing, we need to do is create a contrasts vector called `contrast_kd` for the Mov10_knockdown comparison to control.

```r
## Define contrasts
contrast_kd <- list( "sampletypeMOV10_knockdown", "sampletypecontrol")
```

Use that contrasts vector to extract a results table and store that to a variable called `res_tableKD`.  

```r
# Extract results table
res_tableKD <- results(dds, contrast=contrast_kd)
```

Take a quick peek at the results table containing Wald test statistics for the Control-Knockdown comparison we are interested in and make sure that format is similar to what we observed with the OE.

***

**Exercise**

Create a contrasts vector for the Mov10_overexpression comparison to *all other samples*.

***


> **NOTE: on p-values set to NA**
> > 
> 1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
> 2. If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance. 
> 3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA. 
