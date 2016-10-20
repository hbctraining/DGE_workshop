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


## Visualizing the results

When we are working with large amounts of data it can be useful to display that information graphically to gain more insight. Visualization deserves an entire course of its own (there is that much to know!). During this lesson we will get you started with some basic plots commonly used when exploring differntial gene expression data.

One way to visualize results would be to simply plot the expression data for a handful of our top genes. We could do that by picking out specific genes of interest, for example Mov10:

	# Plot expression for single gene
	plotCounts(dds, gene="MOV10", intgroup="sampletype")
	
![topgene](../img/topgen_plot.png)

### Volcano plot

This would be great to validate a few genes, but for more of a global view there are other plots we can draw. A commonly used one is a volcano plot; in which you have the log transformed adjusted p-values plotted on the y-axis and log2 fold change values on the x-axis. There is no built-in function for the volcano plot in DESeq2, but we can easily draw it using `ggplot2`. First, we will need to create a `data.frame` object from the results, which is currently stored in a `DESeqResults`  object:

	# Create dataframe for plotting
	df <- data.frame(res_tableOE)

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

Alternatively, we could extract only the genes that are identifed as significant and the plot the expression of those genes using a heatmap.


First, let's sort the results file by adjusted p-value:
	
	### Sort the results tables
	res_tableOE_sorted <- res_tableOE[order(res_tableOE$padj), ]
	res_tableKD_sorted <- res_tableKD[order(res_tableKD$padj), ]
	
Now let's get the gene names for those significant genes:

	### Get significant genes
	sigOE <- row.names(res_tableOE_sorted)[which(res_tableOE_sorted$threshold)]
	sigKD <- row.names(res_tableKD_sorted)[which(res_tableKD_sorted$threshold)]
	
We can then use those genes to select the corresponding rows from the normalized data matrix:

	### Extract normalized expression for significant genes
	norm_OEsig <- normalized_counts[sigOE,]

Now let's draw the heatmap using `pheatmap`:

	### Annotate our heatmap (optional)
	annotation <- data.frame(sampletype=meta[,'sampletype'], 
                         row.names=row.names(meta))

	### Set a color palette
	heat.colors <- brewer.pal(6, "YlOrRd")
	
	### Run pheatmap
	pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,
	annotation= annotation, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)
         
![sigOE_heatmap](../img/sigOE_heatmap.png)       

> *NOTE:* There are several additional arguments we have included in the function for aesthetics. One important one is `scale="row"`, in which Z-scores are plotted, rather than the actual normalized count value. Z-scores are computed on a gene-by-gene basis by subtracting the mean and then dividing by the standard deviation. The Z-scores are computed **after the clustering**, so that it only affects the graphical aesthetics and the color visualization is improved.

***

**Exercise**

1. Generate two figures for the KD-control comparison: a volcano plot and a heatmap. 
2. Save both images to file.

***



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*
