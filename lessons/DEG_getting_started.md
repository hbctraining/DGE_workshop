---
title: "Setting up for Gene-level differential expression analysis"
author: "Meeta Mistry, Radhika Khetani"
date: "April 26, 2017"
---

Approximate time: 40 minutes

## Learning Objectives 

* Have a general idea of the experiment and its objectives
* Understand how and why we choose this dataset
* Getting setup in R (project setup, loading data, loading libraries)

 
## Understanding the dataset

We will be using a real RNA-Seq dataset for today's class. It is part of a larger study described in [Kenny PJ et al, Cell Rep 2014](http://www.ncbi.nlm.nih.gov/pubmed/25464849). 

We are only using the [RNA-Seq](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499) dataset which is publicly available in the [SRA](http://www.ncbi.nlm.nih.gov/sra). The RNA-Seq was performed on HEK293F cells that were either transfected with a MOV10 transgene, or siRNA to knock down Mov10 expression, or non-specific (irrelevant) siRNA. This resulted in 3 conditions **Mov10 oe** (over expression), **Mov10 kd** (knock down) and **Irrelevant kd**, respectively. The number of replicates is as shown below. 

Using these data, we will evaluate transcriptional patterns associated with perturbation of MOV10 expression. Please note that the irrelevant siRNA will be treated as our control condition.

<img src="../img/dataset.png" width="400">



***What is the purpose of these datasets? What does Mov10 do?***

The authors are investigating interactions between various genes involved in Fragile X syndrome, a disease in which there is aberrant production of the FMRP protein. 

> **FMRP** is “most commonly found in the brain, is essential for normal cognitive development and female reproductive function. Mutations of this gene can lead to fragile X syndrome, mental retardation, premature ovarian failure, autism, Parkinson's disease, developmental delays and other cognitive deficits.” - from [wikipedia](https://en.wikipedia.org/wiki/FMR1)

> **MOV10**, is a putative RNA helicase that is also associated with **FMRP** in the context of the microRNA pathway. 

**The hypothesis [the paper](http://www.ncbi.nlm.nih.gov/pubmed/25464849) is testing is that FMRP and MOV10 associate and regulate the translation of a subset of RNAs.**

<img src="../img/mov10-model.png" width="400">

**Our questions:**
* What patterns of expression can we identify with the loss or gain of MOV10? 
* Are there any genes shared between the two conditions?

## Metadata

In addition to the raw sequence data that is available in SRA we also need to collect **information about the data**, also known as **metadata**.

Data sharing is important in the biological sciences to promote scientific integrity, and disseminate scientific discovery; but it can be difficult if all of the required information is not provided. From the SRA we can retrieve the sequence data (FASTQ files), but how useful is it if we know nothing about the samples that this sequence data originated from? **Metadata is a broadly used term which encompasses any kind of information that relates to our data, whether it is about the experimental design (i.e genotype) or metrics related to the sequence data (i.e sequencing depth).**

Below is some of the metadata associated with the dataset we are using today.

* The RNA was extracted from treated **HEK293F cells**.  
* The cDNA libraries for this dataset are **stranded** and were generated using the **TruSeq Stranded mRNA Library Prep Kit** from Illumian. 
* Sequencing was carried out on the **Illumina HiSeq-2500 for 100bp single end** reads. 
* **~40 million reads** per sample were generated.

> [Metadata generation/maintenance](http://datamanagement.hms.harvard.edu/metadata-overview) is part of ["Data Management"](http://datamanagement.hms.harvard.edu/biomedical-data-management-planning), which is an important aspect of working with large datasets.

***

**Exercise**

1. What types of metadata are used in your experimental design (any experiment)?
2. What other kinds of metadata might a sequencing project generate?
3. Why is this type of information important?

***


## Setting up

Let's get started by opening up RStudio and setting up a new project for this analysis. 

1. Go to the `File` menu and select `New Project`.
2. In the `New Project` window, choose `New Directory`. Then, choose `Empty Project`. Name your new directory `DEanalysis` and then "Create the project as subdirectory of:" the Desktop (or location of your choice).
3. The new project should automatically open in RStudio. 

To check whether or not you are in the correct working directory, use `getwd()`. It shoud return a path similar to `/../../DEanalysis` in the console. Create 3 new folders/directories called `data`, `meta` and `results` using the `New folder` button in the `Files` tab (panel on the right side under the Environment panel). Remember the key to a good analysis is keeping organized from the start!

Go to the `File` menu at the top left, and select `New File` followed by `R Script`. This should open up a script editor in the top left hand corner. This is where we will be typing and saving all commands required for this analysis, just as we did the R course. Let's get started by typing in the following (commented) header lines:

```
## Gene-level differential expression analysis using DESeq2
## May 18th, 2017
```

Now save the file as `de_script.R`. When finished your working directory should now look similar to this:

![setup](../img/settingup.png)

Finally, we need to grab the files that we will be working with for the analysis. Right click on the links below, and choose the "Save link as ..." option to download:

* Save the [full counts matrix](https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_counts.txt) file in the `data` directory.
* Save the [full metadata table](https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_meta.txt) file in the `meta` directory.

## Loading libraries

For this analysis we will be using several R packages, some which have been installed from CRAN and others from Bioconductor. To use these packages (and the functions contained within them), we need to **load the libraries.** Add the following to your script and don't forget to comment liberally!

```r
## Setup
### Bioconductor and CRAN libraries used
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```

## Loading data

To load the data into our current environment, we will be using the `read.table` function. We need to provide the path to each file and also specify arguments to let R know that we have a header (`header = T`) and the first column is our row names (`row.names =1`). By default the function expects tab-delimited files, which is what we have.

```r
## Load in data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)
```

Use `class()` to inspect our data and make sure we are working with data frames:

```r
### Check classes of the data we just brought in
class(data)
class(meta)
```

## View your data

Make sure your datasets contain the expected samples / information before proceeding to perfom any type of analysis. 

```r
View(data)
View(meta)
```

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
