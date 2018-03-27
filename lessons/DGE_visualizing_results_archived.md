---
title: "DEG visualization of results"
author: "Meeta Mistry"
date: "October 17, 2016"
---

Approximate time: 75 minutes

## Learning Objectives 

* Exploring our significant genes using data visualization
* Using volcano plots to evaluate relationships between DEG statistics
* Plotting expression of significant genes using heatmaps

## Extracting significant genes

The FDR threshold on it's own doesn't appear to be reducing the number of significant genes. With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also **add a fold change threshold**. The `summary()` function doesn't have an argument for fold change threshold,

> *NOTE:* the `results()` function does have an option to add a fold change threshold and subset the data this way. Take a look at the help manual using `?results` and see what argument would be required. However, rather than subsetting the results, we want to return the whole dataset and simply identify which genes meet our criteria. 

Let's first create variables that contain our threshold criteria:

```r
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
```

The `lfc.cutoff` is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable. Let's create vector that helps us identify the genes that meet our criteria:

```r
threshold <- res_tableOE$padj < padj.cutoff & abs(res_tableOE$log2FoldChange) > lfc.cutoff
```

We now have a logical vector of values that has a length which is equal to the total number of genes in the dataset. The elements that have a `TRUE` value correspond to genes that meet the criteria (and `FALSE` means it fails). **How many genes are differentially expressed in the Overexpression compared to Control, given our criteria specified above?** Does this reduce our results? 

```r
length(which(threshold))
```
	
To add this vector to our results table we can use the `$` notation to create the column on the left hand side of the assignment operator, and the assign the vector to it instead of using `cbind()`:

```r
res_tableOE$threshold <- threshold                
```

Now we can easily subset the results table to only include those that are significant using the `subset()` function:

```r
sigOE <- data.frame(subset(res_tableOE, threshold == TRUE))
```

Using the same thresholds as above (`padj.cutoff < 0.05` and `lfc.cutoff = 0.58`), create a threshold vector to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.

```r
threshold_KD <- res_tableKD$padj < padj.cutoff & abs(res_tableKD$log2FoldChange) > lfc.cutoff
```

Take this new threshold vector and add it as a new column called `threshold` to the `res_tableKD` which contains a logical vector denoting genes as being differentially expressed or not. **How many genes are differentially expressed in the Knockdown compared to Control?** Subset the data to keep only the significant genes.

```r
res_tableKD$threshold <- threshold_KD

sigKD <- data.frame(subset(res_tableKD, threshold == TRUE))
``` 


## Visualizing the results

When we are working with large amounts of data it can be useful to display that information graphically to gain more insight. Visualization deserves an entire course of its own, but during this lesson we will get you started with some basic plots commonly used when exploring differential gene expression data.

One way to visualize results would be to simply plot the expression data for a handful of our top genes. We could do that by picking out specific genes of interest, for example Mov10:

```r
# Plot expression for single gene
plotCounts(dds, gene="MOV10", intgroup="sampletype")
```
	
![topgene](../img/topgen_plot.png)

### Volcano plot

The above plot would be great to validate a select few genes, but for more of a global view there are other plots we can draw. A commonly used one is a volcano plot; in which you have the log transformed adjusted p-values plotted on the y-axis and log2 fold change values on the x-axis. There is no built-in function for the volcano plot in DESeq2, but we can easily draw it using `ggplot2`. First, we will need to create a `data.frame` object from the results, which is currently stored in a `DESeqResults`  object:

```r
# Create dataframe for plotting
df <- data.frame(res_tableOE)
```

Now we can start plotting. The `geom_point` object is most applicable, as this is essentially a scatter plot:

```
# Volcano plot
ggplot(df) +
	geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
	xlim(c(-2,2)) +
	ggtitle('Mov10 overexpression') +
	xlab("log2 fold change") + 
	ylab("-log10 adjusted p-value") +
	theme(legend.position = "none",
    	plot.title = element_text(size = rel(1.5)),
    	axis.title = element_text(size = rel(1.5)),
    	axis.text = element_text(size = rel(1.25)))  
```

![volcano](../img/volcanoplot-1.png)


### Heatmap

Alternatively, we could extract only the genes that are identified as significant and the plot the expression of those genes using a heatmap. We can then use the genes from the subsetted data frames to select the corresponding rows from the normalized data matrix:

```r
### Extract normalized expression for significant genes
norm_OEsig <- normalized_counts[rownames(sigOE),]
```

Now let's draw the heatmap using `pheatmap`:

```r
### Annotate our heatmap (optional)
annotation <- data.frame(sampletype=meta[,'sampletype'], 
                     row.names=row.names(meta))

### Set a color palette
heat.colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,
annotation= annotation, border_color=NA, fontsize = 10, scale="row",
     fontsize_row = 10, height=20)
```
         
![sigOE_heatmap](../img/sigOE_heatmap.png)       

> *NOTE:* There are several additional arguments we have included in the function for aesthetics. One important one is `scale="row"`, in which Z-scores are plotted, rather than the actual normalized count value. Z-scores are computed on a gene-by-gene basis by subtracting the mean and then dividing by the standard deviation. The Z-scores are computed **after the clustering**, so that it only affects the graphical aesthetics and the color visualization is improved.

***

**Exercise**

1. Generate two figures for the KD-control comparison: a volcano plot and a heatmap. 
2. Save both images to file.

***
->**NOTE:** Advanced visualization methods for DE analysis results are available as [additional material](https://github.com/hbctraining/Training-modules/blob/master/Visualization_in_R/lessons/03_advanced_visualizations.md).

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*
