---
title: "Quantitation of transcript abundance using Salmon"
author: "Mary Piper and Meeta Mistry"
date: "Wednesday, September 20, 2017"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 30 minutes

## Learning Objectives

* Explore using pseudocounts from lightweight algorithms to perform differential gene expression analysis
* Understand different DESeq2 functions when working with pseudocounts


## Lightweight alignment and quantification of gene expression

In the standard RNA-seq pipeline, we obtain a matrix of raw counts which represent the total number of reads that map to a gene in a particular sample. **Another strategy for quantification which has more recently been introduced involves transcriptome mapping**. Tools that fall in this category include [Kallisto](https://pachterlab.github.io/kallisto/about), [Sailfish](http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html) and [Salmon](https://combine-lab.github.io/salmon/); each working slightly different from one another. (For this workshop we will explore Salmon pseudocounts in more detail.) Common to all of these tools is that **base-to-base alignment of the reads is avoided**, which is a time-consuming step, and these tools **provide quantification estimates much faster than do standard approaches** (typically more than 20 times faster) with **improvements in accuracy** at **the transcript level**. 

These transcript expression estimates, often referred to as **'pseudocounts'**, can be converted for use with DGE tools like DESeq2 or the estimates can be used directly for isoform-level differential expression using a tool like [Sleuth](http://www.biorxiv.org/content/biorxiv/early/2016/06/10/058164.full.pdf). 

<img src="../img/alignmentfree_workflow2.png" width="500">


### What is Salmon?

[Salmon](http://salmon.readthedocs.io/en/latest/salmon.html#using-salmon) is based on the philosophy of lightweight algorithms, which use the reference transcriptome (in FASTA format) and raw sequencing reads (in FASTQ format) as input, but do not align the full reads. These tools perform both mapping and quantification. Unlike most lightweight and standard alignment/quantification tools, **Salmon utilizes sample-specific bias models for transcriptome-wide abundance estimation**. Sample-specific bias models are helpful when needing to account for known biases present in RNA-Seq data including:

- GC bias
- positional coverage biases
- sequence biases at 5' and 3' ends of the fragments
- fragment length distribution
- strand-specific methods

If not accounted for, these biases can lead to unacceptable false positive rates in differential expression studies [[2](http://salmon.readthedocs.io/en/latest/salmon.html#quasi-mapping-based-mode-including-lightweight-alignment)]. The **Salmon algorithm can learn these sample-specific biases and account for them in the transcript abundance estimates**. Salmon is extremely fast at "mapping" reads to the transcriptome and often more accurate than standard approaches [[2](http://salmon.readthedocs.io/en/latest/salmon.html#quasi-mapping-based-mode-including-lightweight-alignment)]. 



### Salmon output

You should see a new directory has been created that is named by the string value you provided in the `-o` command. Take a look at what is contained in this directory:

    $ ls -l Mov10_oe_1.subset.salmon/
    
There is a logs directory, which contains all of the text that was printed to screen as Sailfish was running. Additionally, there is a file called `quant.sf`. 

This is the **quantification file** in which each row corresponds to a transcript, listed by Ensembl ID, and the columns correspond to metrics for each transcript:

```bash
Name    Length  EffectiveLength TPM     NumReads
ENST00000632684.1       12      3.00168 0       0
ENST00000434970.2       9       2.81792 0       0
ENST00000448914.1       13      3.04008 0       0
ENST00000415118.1       8       2.72193 0       0
ENST00000631435.1       12      3.00168 0       0
ENST00000390567.1       20      3.18453 0       0
ENST00000439842.1       11      2.95387 0       0

....

```

*  The first two columns are self-explanatory, the **name** of the transcript and the **length of the transcript** in base pairs (bp). 
*  The **effective length** represents the the various factors that effect the length of transcript due to technical limitations of the sequencing platform.
* Salmon outputs ‘pseudocounts’ which predict the relative abundance of different isoforms in the form of three possible metrics (KPKM, RPKM, and TPM). **TPM (transcripts per million)** is a commonly used normalization method as described in [[1]](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2820677/) and is computed based on the effective length of the transcript.
* Estimated **number of reads** (an estimate of the number of reads drawn from this transcript given the transcript’s relative abundance and length)


## Next steps:  Performing DE analysis on Pseudocounts

The pseudocounts generated by Salmon are represented as normalized TPM (transcripts per million) counts and map to transcripts or genes, depending on the input in the index step; in our case it was transcripts. These **need to be converted into non-normalized count estimates for performing DE analysis using existing count-based tools**. We will not cover this in the workshop but have [additional material](DE_analysis.md) to walk you through it if you are interested. 


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*





