---
title: "Differential expression of transcripts using Sleuth"
author: "Mary Piper"
date: "Tuesday, March 27th, 2018"
---

Contributors: Mary Piper

Approximate time: 1.25 hours

## Learning Objectives

* Determine differential expression of isoforms and visualization of results using Sleuth
* Understand how Sleuth determines biological and technical variation

## Transcript-level differential expression

Until this point we have focused on looking for expression changes at the gene-level. However, if you are interested in looking at **splice isoform expression changes between groups** the previous methods (i.e DESeq2) will not work. To demonstrate how to identify transcript-level differential expression we will be using a tool called Sleuth.

<img src="../img/sleuth_fullworkflow_new.png" width="500">

## What is Sleuth?

[Sleuth](http://pachterlab.github.io/sleuth/) is a fast, lightweight tool that uses transcript abundance estimates output from pseudo-alignment algorithms that use bootstrap sampling, such as Sailfish, Salmon, and Kallisto, to perform differential expression analysis of gene isoforms. 

To analyze the differential expression of gene isoforms, it is expected that RNA-Seq reads will often align to multiple isoforms of the same gene. Therefore, **multimapping reads cannot be ignored** to properly determine abundances of gene isoforms. 

Due to the statistical procedure required to assign reads to gene isoforms, in addition to the random processes underlying RNA-Seq, there will be **technical variability in the abundance estimates** output from the pseudo-alignment tool [[2](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html), [3](http://biorxiv.org/content/biorxiv/early/2016/06/10/058164.full.pdf)]. For example, if we performed multiple technical replicates and estimated abundances for gene isoforms, the abundance estimates for the technical replicates would exhibit variability greater than expected due to variability in the transcript abundance estimation process. Therefore, **we would need technical replicates to distinguish technical variability from the biological variability**.

Sleuth accounts for this technical variability by using **bootstraps as a proxy for technical replicates**, which are used to model the variability in the abundance estimates. Bootstrapping essentially **calculates the abundance estimates for all genes using a different sub-sample of reads** during each round of bootstrapping. The variation in the abundance estimates output from each round of bootstrapping is used for the estimation of the technical variance for each gene. 

<p align="center">
  <img src="../img/sleuth_tech_var.png" width="600"/>
</p>

_Adapted from: Nature Methods **14**, 687–690 (2017)_

Sleuth models the unobserved true abundance (logarithm of true counts) using a general linear model, but includes the technical variance (variance between bootstrapping runs) as error in the response variable. 

<p align="center">
  <img src="../img/sleuth_formula1.png" width="500"/>
</p>


The observed (log) abundance estimates represent the sum of the true counts and the technical noise. Thereby, sleuth teases apart the source of the variance (technical vs. biological) for estimation of the "true" biological variance when determining whether transcripts are differentially expressed.

<p align="center">
  <img src="../img/sleuth_formula2.png" width="425"/>
</p>

In addition to performing differential expression analysis of transcripts, the sleuth tool also provides an html interface allowing exploration of the data and differential expression results interactively. More information about the theory/process for sleuth is available in the [Nature Methods paper](https://www-nature-com.ezp-prod1.hul.harvard.edu/articles/nmeth.4324), this [blogpost](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/) and step-by-step tutorials are available on the [sleuth website](https://pachterlab.github.io/sleuth/walkthroughs).

***NOTE:*** *Kallisto is distributed under a non-commercial license, while Sailfish and Salmon are distributed under the [GNU General Public License, version 3](http://www.gnu.org/licenses/gpl.html).*


## Set-up for Running Sleuth

Sleuth is a lightweight algorithm that can be quickly run on our personal computers [[2](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html)]. We can use the previous `DE_pseudocounts` project we created previously.

1. Open the RStudio project entitled `DE_pseudocounts`.
2. Create a new R script ('File' -> 'New File' -> 'Rscript'), and save it as `sleuth_de.R`

To perform any analysis, we need to load the libraries for `wasabi` and `sleuth`. Sleuth also requires `annotables`, so this package will be loaded as well:

```R
library(wasabi)
library(sleuth)
library(annotables)
library(tidyverse)
```

## Using Wasabi to convert Salmon output for Sleuth

Sleuth was built to use the bootstrapped estimates of transcript abundance from Kallisto; however, abundance estimates from Salmon (or Sailfish) work just as well, as long as bootstrapping is performed. To use our Salmon output estimates, we need to convert them to a Sleuth-compatible format using the Wasabi package.

First, we create a simple vector containing the paths to the directories containing the transcript abundance estimates for each sample (folders containing the .quant files). We can use the `file.path()` function to give the paths to each of the directories. 

Now, let's use this function to create our list of the paths to our transcript abundance files:

```R
sf_dirs <- file.path("data", "salmon", c("Mov10_kd_2.salmon", "Mov10_kd_3.salmon", "Mov10_oe_1.salmon", "Mov10_oe_2.salmon", "Mov10_oe_3.salmon","Irrel_kd_1.salmon", "Irrel_kd_2.salmon", "Irrel_kd_3.salmon"))

sf_dirs
```

Now, we simply run the `prepare_fish_for_sleuth` function, which will write some status messages to the console and, when it's done, each directory will contain an `abundance.h5` file in a sleuth-compatible format.

```R
prepare_fish_for_sleuth(sf_dirs)
```

Each of the sample directories should now contain the `abundance.h5` files. These 
files will be used as input to Sleuth.

<p align="center">
  <img src="../img/sleuth_abund_h5.png" width="600"/>
</p>

## Sleuth for estimation of differential expression of transcripts

![sleuth](../img/sleuth_workflow.png)

The workflow for Sleuth is similar to the workflow followed for DESeq2, even though, the models for estimating differential expression are very different. 

**Step 1:** Creation of Sleuth object to provide metadata, estimated counts, and design formula for the analysis, in addition to a biomaRt database to switch between transcript IDs and associated gene names.

**Step 2:** Fit the sleuth model
	
- Estimation of size (normalization) factors using the median of ratios method (similar to DESeq2)
	
- Normalization of estimated counts using size factors (est. counts / size factors - similar to DESeq2)

- Filtering of low abundance transcripts (< 5 est counts in more than 47% of the samples)

- Normalization of technical variation estimates

- Estimation of biological variance and shrinkage estimates (With small sample sizes, we will make very bad estimates of transcript-wise dispersion unless we share information across transcripts. Sleuth regularizes the biological variance estimate with shrinkage, similar to DESeq2, except uses a different statistical method (similar to Limma Voom).)
	
- Parameter estimation and estimation of variance using the general linear model.

- Identification of:
	- **Coefficients:** indicating overall expression strength
	- **Beta values:** estimating fold changes	

**Step 3:** Test for significant differences between conditions

After performing all analysis steps, we will explore our results by transferring the results to our local machine, and we will use the html interface available through the sleuth package.

## Sleuth workflow

### Step 1: Create Sleuth object for analysis

![sleuth](../img/sleuth_workflow1.png)

Similar to DESeq2, we need to tell Sleuth where to find the **metadata** (specifying which samplegroups the samples belong to, and any other metadata we want included in the analysis), **estimated counts** (output from Salmon) and the **design formula**. In addition, we also need a **biomaRt database** to easily convert between transcript IDs and associated gene names. To create this object there is no simple function like in DESeq2 (e.g. DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)). 

To create this Sleuth object, we need to perform the following steps:

1. Create a dataframe containing metadata and locations of the estimated counts files:

  	- including any columns containing metadata to used in the analysis
  	- a column named `sample` containing all of the sample names matching the names in the metadata file
  	- a column named `path` containing the path to the abundance estimate files output from `wasabi`
        
2. Create a variable containing the model design 

3. Use biomaRt to create a dataset for Sleuth to query for Ensembl IDs and associated gene names

#### Create a dataframe needed to generate Sleuth analysis object

Read in the metadata file and use the `data.frame()` function to ensure it is a dataframe, then combine the metadata with the paths to the transcript abundance files to use as input for the Sleuth analysis. 

First we need to read in the metadata file:

```r
# Read in metadata file

summarydata <- data.frame(read.table("meta/Mov10_full_meta.txt", header=TRUE, row.names=1), check.rows=FALSE)

summarydata
```

Then we make sure the metadata and count estimate sample names match:

```r
# Make sure the order of the `sfdirs` created above matches the order of samples in the `summarydata` rownames

sf_dirs_samples <- sf_dirs %>%
  basename %>% 
  str_replace(pattern = "\\.salmon", "")

all(sf_dirs_samples == rownames(summarydata))
```
Now, we can name the vector of directory paths with the corresponding sample names

```r
# Name the directory paths for the abundance files with their corresponding sample IDs

names(sf_dirs) <- rownames(summarydata)

sf_dirs
```

Finally, we can generate the data frame containing the metadata:

```r
# Generate the dataframe

sfdata <- summarydata
```

Sleuth expects the data to be presented in a specific format with specific column and row names; therefore, we will need to name columns based on the sleuth requirements for the analysis. 

Sleuth requires a column entitled "sample" containing the sample names:

```r
# Adding a column named 'sample'

sfdata$sample <- rownames(sfdata)
```

Now, we can include the path to the count estimate folders. Sleuth requires a column entitled "path" containing the paths to the estimated counts files stored in our `sf_dirs`:

```r
sfdata$path <- sf_dirs

sfdata
```

#### Provide the model design

Now that we have the metadata and location of the count estimates, we can input our design formula to determine the covariates and/or confounders that should be included in your experimental design model. Sleuth can be used to analyze multiple conditions from complex experimental designs.

Within Sleuth, models are written similar to DESeq2. Since the only condition we plan to test is our sample type, our design formula is very simple:

```r
design <- ~ sampletype
```

More complex designs can be analyzed using Sleuth as well by adding additional covariates (i.e `design <- ~ sex + treatment`). Interaction terms can also be added to the design formula to test if the effect attributable to a given condition is different based on another factor, for example, if the treatment effect differs between sexes. To learn more about setting up design formulas for more complex designs, see the [DESeq2 tutorial](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) (chapter 3 discusses complex experimental designs). There is also a [recent post](http://nxn.se/post/134227694720/timecourse-analysis-with-sleuth) describing the use of Sleuth to perform time course analyses. While Sleuth has much flexiblity in design models, it is unable to support some complex designs, such as nested models.

#### Create Biomart dataset to query

The last component to include for our analysis is the biomaRt Ensembl genome database to obtain the Ensembl transcript/gene IDs and gene names for annotation of results. BiomaRt allows extensive genome information to be accessible during an analysis.

```r
# Using annotables
t2g <- grch37_tx2gene

t2g <- merge(x= grch37[, c("symbol", "ensgene")], y = t2g, by.x="ensgene", by.y= "ensgene")

## Rename the columns for use in Sleuth

t2g <- dplyr::rename(t2g, target_id = enstxp , 
                     ens_gene = ensgene, 
                     ext_gene = symbol)
```

### Step 2: Fit the sleuth model

![sleuth](../img/sleuth_workflow2.png)

#### Fit the transcript abundance data to the Sleuth model

```r
# Create sleuth object for analysis 

so <- sleuth_prep(sfdata, design, target_mapping = t2g, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)  

# Fit the transcript abundance data to the sleuth model

so <- sleuth_fit(so)

# NOTE: alternatively the two prior steps could have been run as: "so <- sleuth_prep(sfdata, design, target_mapping = t2g) %>% sleuth_fit()

```

#### Check which models have been fit and which coefficients can be tested

Ensure the design model and coefficients are correct for your analysis. The level not shown is the base level.

```r
models(so)
```
> **NOTE:** Sleuth will automatically use the first level (alphabetically) in the factor variable being tested to compare all other conditions against (in our metadata, this is 'control'). If you want to use a different condition to be the base level, then you would need to use the relevel() function to change the base level of the variable in step 1 above. For example, if we wanted the base level of `sampletype` to be "MOV10_knockdown", we could use the following code:
>
>```r
> # DO NOT RUN!
> summarydata$sampletype <- relevel(summarydata$sampletype, ref = "MOV10_knockdown")
>```
>***An ordered factor will not give interpretable output, so do not order the factor using the factor() function, use relevel() instead.***

### Step 3: Test significant differences between conditions using the Wald test

![sleuth](../img/sleuth_workflow3.png)

At this step in the workflow, we need to specify which level we want to compare against the base level (use the name given for the coefficients from `models(so)`):

```r
# Wald test for differential expression of isoforms

oe <- sleuth_wt(so, 'sampletypeMOV10_overexpression')

# output results

sleuth_results_oe <- sleuth_results(oe, 'sampletypeMOV10_overexpression', show_all = TRUE)
```

>**NOTE:** There are also methods for performing the LRT test and specifying a full and reduced model, which are described in detail in a [sleuth walk-through](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html).

The output represents the results from the differential expression testing.

![sleuth_results](../img/sleuth_mov10_results.png)

### Exploring transcript-level expression between samples

Now we can perform some exploratory analyses, such as PCA and heatmap:

```r
plot_pca(oe, color_by = "sampletype")

plot_sample_heatmap(oe)
```
<p align="center">
  <img src="../img/sleuth_mov10_pca.png" width="400"/>
</p>

<p align="center">
  <img src="../img/sleuth_mov10_heatmap.png" width="400"/>
</p>

In addition there is a histogram plot to explore count distributions between sample groups, which should be similar for DE testing. The count distributions represent the proportion of genes (on the y-axis) associated with the number of counts (designated on the x-axis):

```r
plot_group_density(oe, use_filtered = FALSE, units = "est_counts",
                   trans = "log", grouping = "sampletype")
```

<p align="center">
  <img src="../img/sleuth_mov10_density_no_filter.png" width="400"/>
</p>

As we know, most genes have few counts, but we filter these genes prior to performing DE analysis. If we want to look at the distributions of the filtered genes used for DE analysis, we could change the `use_filtered` argument to `TRUE`.

```r
plot_group_density(oe, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = "sampletype")
```

<p align="center">
  <img src="../img/sleuth_mov10_density.png" width="400"/>
</p>

There are also functions to explore the results, such as the MA plot:

```r
# Try the plot_ma() function
plot_ma(oe) # shows that we need to specify the 'test'

?plot_ma

# Check for possible tests
tests(oe)

# Try the plot_ma() function again
plot_ma(oe, test="sampletypeMOV10_overexpression", sig_level = 0.05)
```

<p align="center">
  <img src="../img/sleuth_mov10_ma.png" width="400"/>
</p>

Sleuth also has some handy functions to plot expression of transcripts with bootstrap variation:

```r
# Plotting

plot_bootstrap(oe, "ENST00000495374.5", units = "est_counts", color_by = "sampletype")

plot_bootstrap(oe, "ENST00000367412.1", units = "est_counts", color_by = "sampletype")
```

<p align="center">
  <img src="../img/mov10_bootstraps.png" width="400"/>
</p>

<p align="center">
  <img src="../img/sleuth_mov10_bootstraps2.png" width="400"/>
</p>

Sleuth also offers us the option to explore the data and results interactively using a web interface. 

```r
sleuth_live(oe)
```

Using this web interface, we can explore diagnostic plots and summaries of experimental factors and expression levels. We also have the ability to perform clustering analyses such as PCA and heatmaps. Finally, we can analyze the differential expression results by plotting MA and volcano plots and by exploring expression levels at the transcript and gene levels. A brief tutorial explaining the options available interactively can be found [here](http://www.rna-seqblog.com/iihg-intro-to-sleuth-for-rna-seq/).

Explore the results table under the `analyses` tab, `test table`.

![results_table](../img/mov10_test_table.png)

Look at the expression levels of Mov10 for three different isoforms using the `analyses` tab, `gene view`. Look for the Mov10 gene, ENSG00000155363. Compare expression of the different isoforms of Mov10.

![mov_isoforms](../img/mov10_isoform_expression.png)

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
