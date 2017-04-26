# Differential gene expression analysis overview 

The goal of RNA-Seq is often to perform differential expression testing to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the condition(s) of interest. 

To determine the expression levels of genes, an RNA-Seq workflow is followed with the steps detailed in the image below. All steps are performed on the command line (Linux/Unix) through the generation of the read counts per gene. The differential expression analysis and any downstream functional analysis are generally performed in R using R packages specifically designed for the complex statistical analyses required to determine whether genes are differentially expressed.


<img src="../img/rnaseq_full_workflow.png" width="400">

An in-depth explanation of these steps is outside the scope of today's class, but a couple of points:

- Even though this flow diagram only shows 1 tool per step after the sequencing step, there are several tools available.
- This is the more standard workflow with an alignment + a counting step, but more recently people are moving to an alignment-free counting workflow using tools like [Salmon](https://combine-lab.github.io/salmon/) and [Kallisto](https://pachterlab.github.io/kallisto/about.html). These newer tools will generate an abundance estimate for the genes, instead of "raw" counts, but the downstream steps for statistical analysis are similar. 

In the next few lessons, we will walk you through an **end-to-end gene-level RNA-seq differential expression workflow** using various R packages. We will start with the count matrix, perform exploratory data analysis for quality assessment and to explore the relationship between samples, perform differential expression analysis, and visually explore the results.

The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher level of expression of that gene in the sample. 

<img src="../img/deseq_counts_overview.png" width="600">

With differential expression analysis, we are looking for genes that change in expression between two or more groups
- case vs. control
- correlation of expression with some variable or clinical outcome

**Why does it not work to identify differentially expressed gene by ranking the genes by how different they are between the two groups (based on fold change values)?**


<img src="../img/foldchange_heatmap.png" width="200">

There are many sources of variation present in your data, and even though the mean expression levels between groups may be quite different, the variation in the data may be so great that the difference in means are not actually significant. 

<img src="../img/de_variation.png" width="500">

This is illustrated for 'GeneA' expression between 'untreated' and 'treated' groups in the figure below. The mean expression level of geneA for the 'treated' group is twice as large as for the 'untreated' group, but the variation between replicates indicates that this may not be a significant difference. **We need to take the variation in the data into account when determining whether genes are differentially expressed.** 

<img src="../img/de_norm_counts_var.png" width="400">


**Essentially, the goal of differential expression analysis is to determine whether the differences in expression (counts) between groups is significant given the variation within groups (replicates) for each gene.** To test for significance, we need an appropriate statistical model that accurately performs normalization (to account for differences in sequencing depth, etc.) and variance modeling (to account for few numbers of replicates and large dynamic expression range). 



## RNA-Seq count distribution

To determine the appropriate statistical model, we need information about the distribution of counts. To get an idea about how RNA-Seq counts are distributed, let's plot the counts for a single sample, 'Mov10_oe_1':

```r
ggplot(data) +
        geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200)
```

<img src="../img/deseq_counts_distribution.png" width="400">

If we zoom in close to zero, we can see a large number of genes with counts of zero:

```r
ggplot(data) +
        geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) + 
        xlim(-5, 500) 
```

<img src="../img/deseq_counts_distribution_zoomed.png" width="400">

These images illustrate some common features of RNA-Seq count data, including a **low number of counts associated with the a large proportion of genes**, and a long right tail due to the **lack of any upper limit for expression**. Unlike microarray data, which has a dynamic range maximum limited due to when the probes max out, there is no limit of maximum expression for RNA-Seq data. Due to the differences in these technologies, the statistical models used to fit the data are different between the two methods. 

> **NOTE:** The log intensities of the microarray data approximate a normal distribution. However, due to the different properties of the of RNA-Seq count data, such as integer counts instead of continuous measurements and non-normally distributed data, the normal distribution does not accurately model RNA-Seq counts [[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3541212/)].

## Modeling count data

Count data (discrete) is often modeled using the **binomial distribution**, which can give you the **probability of getting a number of heads upon tossing a coin a number of times**. However, not all count data can be fit with the binomial distribution. 

With some events, like the lottery, when **the number of cases is very large (people who buy tickets), but the probability of an event is very small (probability of winning)**, the **Poisson distribution** is used to model these types of count data. In these cases, the number of events (people who win) generally range between 1 and 10. [Details provided by Rafael Irizarry in the EdX class.](https://youtu.be/fxtB8c3u6l8)

**During RNA-Seq, a very large number of RNAs are present and the probability of pulling out a particular transcript is very small.** However, after taking a large sample, the sum of all counts for that transcript is often between 1 and 10. If this applies, then every gene would follow the Poisson distribution. [A nice description of this concept is presented by Rafael Irizarry in the EdX class.](https://youtu.be/HK7WKsL3c2w)

**If the proportions of mRNA stayed exactly constant between biological replicates** for RNA-Seq data, we could expect Poisson distribution. But realistically, biological variation across biological replicates is expected, and, for RNA-Seq data, genes with larger average expression levels tend to have larger observed variances across replicates. This phenomena of 'having different scatter' is known as data *heteroscedasticity*. **The Negative Binomial (NB) model is a good approximation, to account for this extra variability between replicates.**

<img src="../img/deseq_nb.png" width="400">


>**NOTE:** 
>
> - **Biological replicates** represent multiple samples (i.e. RNA from different mice) representing the same sample class
> - **Technical replicates** represent the same sample (i.e. RNA from the same mouse) but with technical steps replicated
> - Usually biological variance is much greater than technical variance, so we do not need to account for technical variance to identify biological differences in expression
> - **Don't spend money on technical replicates - biological replicates are much more useful**
 
**How do I know if my data should be modeled using the Poisson distribution or Negative Binomial distribution?** 

If it's count data, it should fit the negative binomial, as discussed previously. However, it can be helpful to plot the *mean versus the variance* of your data. *Remember for the Poisson model, mean = variance, but for NB, mean < variance.*

Run the following code to plot the mean versus variance for the 'Mov10 overexpression' replicates:

```r
mean_counts <- apply(data[, 1:3], 1, mean)
variance_counts <- apply(data[, 1:3], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
        geom_point(aes(x=mean_counts, y=variance_counts)) + 
        geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
        scale_y_log10() +
        scale_x_log10()
```

<img src="../img/deseq_mean_vs_variance.png" width="600">

Note that in the figure, the variance across replicates tends to be greater than the mean (red line), especially for genes with large mean expression levels. **This is a good indication that our data do not fit the Poisson distribution and we need to account for this increase in variance using the Negative Binomial model (i.e. Poisson will underestimate variability leading to an increase in false positive DE genes).**

### Improving mean estimates (i.e. reducing variance) with biological replicates

The variance or scatter tends to reduce as we increase the number of biological replicates (*variance will approach the Poisson distribution with increasing numbers of replicates*). Standard deviations of averages are smaller than standard deviations of individual observations. **The value of additional replicates is that as you add more data (replicates), you get increasingly precise estimates of group means, and ultimately greater confidence in the ability to distinguish differences between sample classes (i.e. more DE genes).**

The figure below illustrates the relationship between sequencing depth and number of replicates on the number of differentially expressed genes identified [[1](https://academic.oup.com/bioinformatics/article/30/3/301/228651/RNA-seq-differential-expression-studies-more)]. Note that an **increase in the number of replicates tends to return more DE genes than increasing the sequencing depth**. Therefore, generally more replicates are better than higher sequencing depth, with the caveat that higher depth is required for detection of lowly expressed DE genes and for performing isoform-level differential expression. Generally, the minimum sequencing depth recommended is 20-30 million reads per sample, but we have seen good RNA-Seq experiments with 10 million reads if there are a good number of replicates.

<img src="../img/de_replicates_img.png" width="500">


