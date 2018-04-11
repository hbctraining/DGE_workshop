---
title: "DGE Analysis Homework"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "October 24th, 2017"
---

# Using DESeq2 for gene-level differential expression analysis

- The metadata below describes an RNA-seq analysis experiment, in which the metadata table below and associated count matrix have been loaded into R as `meta` and `counts`, respectively. Additionally, all of the appropriate libraries have been loaded for you. Use the information in the table to answer the following questions.  

	**meta**
	
	| |Genotype	|Celltype	|Batch|
	| ------ | ------- | -------- | --- |
	|sample1	|Wt	|typeA	|second |
	|sample2	|Wt	|typeA	|second|
	|sample3	|Wt	|typeA	|first|
	|sample4	|KO	|typeA	|first|
	|sample5	|KO	|typeA	|first|
	|sample6	|KO	|typeA	|second|
	|sample7	|Wt	|typeB	|second|
	|sample8	|Wt	|typeB	|first|
	|sample9	|Wt	|typeB	|second|
	|sample10	|KO	|typeB	|first|
	|sample11	|KO	|typeB	|first|
	|sample12	|KO	|typeB	|second|


	**NOTE: This is an exercise in thinking about running DESeq2. You do not need to run any code in R/RStudio. Refer to the materials/lessons from class to answer the following questions.**

	**a.** Reorder the columns of the `counts` dataset such that `rownames(meta) == colnames(counts)`.
	
	```
	idx <- match(rownames(meta), colnames(counts))
	
	counts <- counts[ , idx]
	```
	
	**b.** Provide the line of code used to create a DESeqDataSet object called `dds` in which `Genotype` is the factor of interest and `Celltype` and `Batch` are other contributing sources of variation in your data.
	
	```r
	dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Celltype + Batch + Genotype)
	```

	**c.** Provide the line of code required to run DESeq2 on `dds`.
	
	```r
	dds <- DESeq(dds)
	```
	
	**d.** Provide the line of code to create a dispersion plot.
	
	```r
	plotDispEsts(dds)
	```
	
	**e.** Provide the line of code to return the results of a Wald test comparison for `Celltype` categories `typeA` versus `typeB` (i.e the fold changes reported should reflect gene expression changes relative to `typeB`).
	 	
	```r
	res_tableCelltype <- results(dds, contrast = c("Celltype", "typeA" , "typeB"))
	```
	
	**f.** Provide the line of code to return the results with log2 fold change shrinkage performed.
	
	```r
	res_tableCelltype <- lfcShrink(dds, contrast=c("Celltype", "typeA" , "typeB"), res=res_tableCelltype)
	```
	
	**g.** Provide the line of code to write the results of the Wald test with shrunken log2 fold changes to file.
	
	```r
	write.table(res_tableCelltype, file="res_Celltype.txt", sep="\t", quote=F, col.names=NA)
	```
	
	**h.** Provide the line of code to subset the results to return those genes with adjusted p-value < 0.05 and logFold2Change > 1.
	
	```r
	filter(res_tableCelltype, padj < 0.05 & abs(log2FoldChange) > 1)
	```
	
# Working with the DESeq2 results table

- Using the de_script.R that we created in class for the differential expression analysis, change the thresholds for adjusted p-value and log2 fold change to the following values:
 
	```r
	padj.cutoff <- 0.01
	
	lfc.cutoff <- 1.5
	```
	
	Using these new cutoffs, perform the following steps:

	**a.** Subset `res_tableOE` to only return those rows that meet the criteria we specified above (adjusted p-values < 0.01 and log fold changes >1.5). Save the subsetted table to a data frame called `sig_table_hw_oe`. Write the code below:
	
	```r
	sig_table_hw_oe <- filter(res_tableOE, padj < 0.01 & abs(log2FoldChange) > 1.5)
	```
	
	**b.** There is a a DESeq2 function that summarizes how many genes are up- and down-regulated using our criteria for `alpha=0.01`. Use this on the `sig_table_hw_oe`. Write the code you would use, and also, list how many genes are up- and down- regulated.
	
	```r
	summary(sig_table_hw_oe)
		out of 13 with nonzero total read count
		adjusted p-value < 0.01
		LFC > 0 (up) : 12, 92% 
		LFC < 0 (down) : 1, 7.7%
	```
	
	**c.** Get the gene names from `sig_table_hw_oe` and save them to a vector called `sigOE_hw`. Write the code below:
	
	```r
	sigOE_hw <- rownames(sig_table_hw_oe)
	```
	
	**d.** Write the `sigOE_hw` vector of gene names to a file called `sigOE_hw.txt` using the `write()` function. Ensure the genes are listed in a single column. Write the code below.
	
	```r
	write(sigOE_hw, file="results/sigOE_hw.txt", ncol=1)
	
	ACHE
	ASCL1
	C1orf95
	H1F0
	HSPA6
	INHBE
	LOC284344
	LPPR3
	MOV10
	MRC2
	NR1D1
	SCRT1
	SNAP25
	```

	 
# Visualizing Results

- For the genes that are differentially expressed in the knockdown versus control comparison (`res_tableKD`), plot an expression heatmap using normalized counts and `pheatmap()` following the instructions below. Write the code you would use to create the heatmap.

	**a.** The heatmap should only include control and knockdown samples. 
	
	```r
	## Extract normalized expression for significant genes
	norm_KDsig <- normalized_counts[rownames(sigKD), c(1,2,6,7,8)]

	## Annotate our heatmap (optional)
	annotation <- meta[c(1,2,6,7,8), "sampletype", drop = FALSE]
	# Or
	annotation <- meta[meta$sampletype != "MOV10_overexpression", "sampletype", drop = FALSE]
	```
	
	**b.** Set up a heat.colors vector using a palette of your choice from brewer.pal (make sure it is different from the one used in class).
	
	```r
	heat.colors <- brewer.pal(6, "YlGnBu")
	```
	
	**c.** Plot the heatmap without clustering the columns. 
	
	```r
	pheatmap(norm_KDsig, color = heat.colors, 
				cluster_rows = T, cluster_cols = F, 
				show_rownames=F, annotation= annotation, 
				fontsize = 10,  fontsize_row = 10, height=20)
	```
	
	**d.** Scale expression values by row.
	
	```r
	pheatmap(norm_KDsig, color = heat.colors, 
				cluster_rows = T, scale="row",
				show_rownames=F, annotation= annotation, 
				fontsize = 10,  fontsize_row = 10, height=20)
	```
	
# Use significant gene lists to find overlaps between the two comparisons 

- Using the original cutoff values, perform the following steps:

	```r
	padj.cutoff <- 0.05

	lfc.cutoff <- 0.58
	```

	**a.** Create separate vectors with gene names for up-regulated genes and down-regulated genes from `res_tableOE` and save as `up_OE` and `down_OE`, respectively. Write the code below:
	
	```r
	up_OE<- row.names(res_tableOE)[which(res_tableOE$padj < padj.cutoff & res_tableOE$log2FoldChange > lfc.cutoff)]

	down_OE <- row.names(res_tableOE)[which(res_tableOE$padj < padj.cutoff & res_tableOE$log2FoldChange < -lfc.cutoff)]
	```
	
	**b.** Create separate vectors with gene names for up-regulated genes and down-regulated genes from `res_tableKD` and save as `up_KD` and `down_KD`, respectively. Write the code below:
	
	```r
	up_KD <- row.names(res_tableKD)[which(res_tableKD$padj < padj.cutoff & res_tableKD$log2FoldChange > lfc.cutoff)]
	
	down_KD <- row.names(res_tableKD)[which(res_tableKD$padj < padj.cutoff & res_tableKD$log2FoldChange < -lfc.cutoff)]
	```
	
	**c.** Test for overlaps between the lists:
	
	- How many, and which genes in `up_OE` are also in `down_KD`?
		
		```r
		length(up_OE[up_OE %in% down_KD])
		# Or
		length(down_KD[down_KD %in% up_OE])
		```
		
		**14 genes**

	- How many, and which genes in `up_KD` are also in `down_OE`?
		
		```r
		length(up_KD[up_KD %in% down_OE])
		# Or
		length(down_OE[down_OE %in% up_KD])
		```
		
		**9 genes**

	
# Using Salmon abundance estimates with DESeq2
	
- In class, we showed materials for using the 'tximport' package in R to summarize the transcript-level estimates generated from Salmon into pseudocounts. 

	Open up RStudio and create a project for using Salmon estimates with DESeq2 (i.e.  ~/Desktop/salmon should be your working directory). 
	
	**a.** Follow the [markdown](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/lessons/DE_analysis.md#differential-expression-analysis-using-pseudocounts) for the Salmon section to create the DESeqDataSet object. 
		
	```r
	# Load libraries
	library(tximport)
	library(readr)
	library(DESeq2)
	library(biomaRt) # tximport requires gene symbols as row names
	
	## List all directories containing data  
	samples <- list.files(path = ".", full.names = F, pattern="\\.salmon$")
    
	## Obtain a vector of all filenames including the path
	files <- file.path(samples, "quant.sf")
    
	## Since all quant files have the same name it is useful to have names for each element
	names(files) <-  samples
	
	tx2gene <- read.delim("tx2gene.txt",sep="\t")
	
	txi <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene, reader=read_tsv, ignoreTxVersion=TRUE)
	
	## Create a sampletable/metadata

	condition=factor(c(rep("Ctl",3), rep("KD", 2), rep("OE", 3)))

	sampleTable <- data.frame(condition, row.names = colnames(txi$counts))

	## Create a DESeqDataSet object
	dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
	```
	
	**b.** Run DESeq2 on the DESeqDataSet object (`dds`) you created.
	
	```r
	dds <- DESeq(dds)
	```
	
	**c.** Use the `results()` function to extract a results table for the OE vs Ctl comparison. Be sure to provide the appropriate contrasts and to perform shrinkage.
	
	```r
	results_OE <- results(dds, contrast = c("condition", "OE", "Ctl"))
	
	results_OE <- lfcShrink(dds, contrast=c("condition", "OE", "Ctl"), res=results_OE)
	```
	
	**d.** Use the `results()` function to extract a results table for the KD vs Ctl comparison. Be sure to provide the appropriate contrasts and to perform shrinkage.
	
	```r
	results_KD <- results(dds, contrast = c("condition", "KD", "Ctl"))
	
	results_KD <- lfcShrink(dds, contrast=c("condition", "KD", "Ctl"), res=results_KD)
	```
	
	**e.** Use the `summary()` function using `alpha=0.05` on the OE results table. 
	
	```r
	summary(results_OE, alpha =  0.05)
	```
	
	- Report the number of genes that are up and down-regulated: 
	
		```r
		LFC > 0 (up)     : 1795, 6.6% 
		
		LFC < 0 (down)   : 2071, 7.6% 
		```

	- Report how many genes are removed due to independent filtering (HINT: these are genes that are filtered due to low mean count):
	
		```r
		low counts [2]   : 11902, 44% 
		```


	- Report how many genes were removed due to an extreme outlier count (HINT: these 'outliers' are identified based on Cook's distance):

	
		```r
		outliers [1]     : 30, 0.11% 
		```


	**f.** Use the `summary()` function using `alpha=0.05` on the KD results table. 

	- Report the number of genes that are up and down-regulated:
	
		```r
		LFC > 0 (up)     : 1554, 5.7%

		LFC < 0 (down)   : 974, 3.6%
		```


	- Report how many genes are removed due to independent filtering (HINT: these are genes that are filtered due to low mean count)
	
		```r
		low counts [2]   : 12940, 48% 
		```


	- Report how many genes were removed due to an extreme outlier count (HINT: these 'outliers' are identified based on Cook's distance)
	
		```r
		outliers [1]     : 30, 0.11%
		```


	**g.** Create a vector called `criteria_oe`, which contains the indexes for which rows in the OE results table have `padj < 0.05` AND `abs(logFC) > 0.58`
	
	```r
	criteria_oe <- which(results_OE$padj < 0.05 & abs(results_OE$log2FoldChange) > 0.58)
	```
	
	**h.** Subset the OE results into a new table using the `criteria_oe` vector:
	
	```r
	results_OE_sig <- results_OE[criteria_oe, ]
	```
	
	**i.**  Create a vector called `criteria_kd`, which contains the indexes for which rows in the KD results table have `padj < 0.05` AND `abs(logFC) > 0.58`
	
	```r
	criteria_kd <- which(results_KD$padj < 0.05 & abs(results_KD$log2FoldChange) > 0.58)
	```
	
	**j.** Subset the KD results into a new table using the `criteria_kd vector`:
	
	```r
	results_KD_sig <- results_KD[criteria_kd, ]
	```
	
	**k.**  Use the following linked datasets as the significant genes using a different workflow (STAR alignment + featureCounts + DESeq2), but same samples. Download them via these links:
	[Download OE genes (870 genes)]
	[Download KD genes (689 genes)]
	Load these gene lists in as `star_KD` and `star_OE`, using the `scan()` function in R (write your code below):
	
	```r
	old_OE <- scan("Mov10_oe_2017.txt", what="character")

	old_KD <- scan("Mov10_kd_2017.txt", what="character")
	```
	
	**l.** Find overlapping OE genes between those that you have identified in your subsetted results table (using the Salmon abundance estimates) to the `star_OE` gene list.
	
	```r
	length(which(row.names(results_OE_sig) %in% old_OE))
	```
	
	- How many genes overlap? **312**

	**m.** Find overlapping KD genes between those that you have identified in your subsetted results table (using the Salmon abundance estimates) to the `star_KD` gene list.
	
	```r
	length(which(row.names(results_KD_sig) %in% old_KD))
	```

	- How many genes overlap? **289**
 
