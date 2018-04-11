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
	
	**b.** Provide the line of code used to create a DESeqDataSet object called `dds` in which `Genotype` is the factor of interest and `Celltype` and `Batch` are other contributing sources of variation in your data.

	**c.** Provide the line of code required to run DESeq2 on `dds`.

	**d.** Provide the line of code to create a dispersion plot.

	**e.** Provide the line of code to return the results of a Wald test comparison for `Celltype` categories `typeA` versus `typeB` (i.e the fold changes reported should reflect gene expression changes relative to `typeB`).
	 
	**f.** Provide the line of code to return the results with log2 fold change shrinkage performed.

	**g.** Provide the line of code to write the results of the Wald test with shrunken log2 fold changes to file.

	**h.** Provide the line of code to subset the results to return those genes with adjusted p-value < 0.05 and logFold2Change > 1.

# Working with the DESeq2 results table

- Using the de_script.R that we created in class for the differential expression analysis, change the thresholds for adjusted p-value and log2 fold change to the following values:
 
	```r
	padj.cutoff <- 0.01
	
	lfc.cutoff <- 1.5
	```
	
	Using these new cutoffs, perform the following steps:

	**a.** Subset `res_tableOE` to only return those rows that meet the criteria we specified above (adjusted p-values < 0.01 and log fold changes >1.5). Save the subsetted table to a data frame called `sig_table_hw_oe`. Write the code below:

	**b.** There is a a DESeq2 function that summarizes how many genes are up- and down-regulated using our criteria for `alpha=0.01`. Use this on the `sig_table_hw_oe`. Write the code you would use, and also, list how many genes are up- and down- regulated.

	**c.** Get the gene names from `sig_table_hw_oe` and save them to a vector called `sigOE_hw`. Write the code below:

	**d.** Write the `sigOE_hw` vector of gene names to a file called `sigOE_hw.txt` using the `write()` function. Ensure the genes are listed in a single column. Write the code below.
	 
# Visualizing Results

- For the genes that are differentially expressed in the knockdown versus control comparison (`res_tableKD`), plot an expression heatmap using normalized counts and `pheatmap()` following the instructions below. Write the code you would use to create the heatmap.

	**a.** The heatmap should only include control and knockdown samples. 

	**b.** Set up a heat.colors vector using a palette of your choice from brewer.pal (make sure it is different from the one used in class).

	**c.** Plot the heatmap without clustering the columns. 

	**d.** Scale expression values by row.

# Use significant gene lists to find overlaps between the two comparisons 

- Using the original cutoff values, perform the following steps:

	```r
	padj.cutoff < 0.05

	lfc.cutoff > 0.58
	```
	
	**a.** Create separate vectors with gene names for up-regulated genes and down-regulated genes from `res_tableOE` and save as `up_OE` and `down_OE`, respectively. Write the code below:

	**b.** Create separate vectors with gene names for up-regulated genes and down-regulated genes from `res_tableKD` and save as `up_KD` and `down_KD`, respectively. Write the code below:

	**c.** Test for overlaps between the lists:
	
	- How many, and which genes in `up_OE` are also in `down_KD`?
	
	- How many, and which genes in `up_KD` are also in `down_OE`?

# Using Salmon abundance estimates with DESeq2
	
- We have materials for using counts from quasialignment tools, such as Salmon, to perform the DE analysis. The 'tximport' package in R is needed to summarize the transcript-level estimates generated from Salmon into pseudocounts. 

	Open up RStudio and create a project for using Salmon estimates with DESeq2 (i.e.  ~/Desktop/salmon should be your working directory). 
	
	**a.** Follow the [markdown](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/lessons/DE_analysis.md#differential-expression-analysis-using-pseudocounts) for the Salmon section to create the DESeqDataSet object. 
	
	**b.** Run DESeq2 on the DESeqDataSet object (`dds`) you created.

	**c.** Use the `results()` function to extract a results table for the OE vs Ctl comparison. Be sure to provide the appropriate contrasts and to perform shrinkage.

	**d.** Use the `results()` function to extract a results table for the KD vs Ctl comparison. Be sure to provide the appropriate contrasts and to perform shrinkage.

	**e.** Use the `summary()` function using `alpha=0.05` on the OE results table. 

	- Report the number of genes that are up and down-regulated:

	- Report how many genes are removed due to independent filtering (HINT: these are genes that are filtered due to low mean count):

	- Report how many genes were removed due to an extreme outlier count (HINT: these 'outliers' are identified based on Cook's distance):


	**f.** Use the `summary()` function using `alpha=0.05` on the KD results table. 

	- Report the number of genes that are up and down-regulated:

	- Report how many genes are removed due to independent filtering (HINT: these are genes that are filtered due to low mean count)

	- Report how many genes were removed due to an extreme outlier count (HINT: these 'outliers' are identified based on Cook's distance)

	**g.** Create a vector called `criteria_oe`, which contains the indexes for which rows in the OE results table have `padj < 0.05` AND `abs(logFC) > 0.58`

	**h.** Subset the OE results into a new table using the `criteria_oe` vector:

	**i.**  Create a vector called `criteria_kd`, which contains the indexes for which rows in the KD results table have `padj < 0.05` AND `abs(logFC) > 0.58`

	**j.** Subset the KD results into a new table using the `criteria_kd vector`:

	**k.**  Use the following linked datasets as the significant genes using a different workflow (STAR alignment + featureCounts + DESeq2), but same samples. Download them via these links:
	
	[Download OE genes (870 genes)](https://wiki.harvard.edu/confluence/download/attachments/216318985/Mov10_oe_2017.txt?version=1&modificationDate=1498507515000&api=v2)
	
	[Download KD genes (689 genes)](https://wiki.harvard.edu/confluence/download/attachments/216318985/Mov10_kd_2017.txt?version=1&modificationDate=1498507515000&api=v2)
	
	Load these gene lists in as `star_KD` and `star_OE`, using the `scan()` function in R (write your code below):

	**l.** Find overlapping OE genes between those that you have identified in your subsetted results table (using the Salmon abundance estimates) to the `star_OE` gene list.
	
	- How many genes overlap?

	**m.** Find overlapping KD genes between those that you have identified in your subsetted results table (using the Salmon abundance estimates) to the `star_KD` gene list.
	
	- How many genes overlap?
 
