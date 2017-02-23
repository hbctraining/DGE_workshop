# Differential gene expression with DESeq2


RNA-Seq Workflow

Goal of differential expression analysis (similar to lecture) -> need appropriate statistical model (normalization and variance modeling). Normalize the counts to make more accurate comparisons between samples of different sequencing depths and variance modeling to account for the increase in variance associated with the few number of samples and large dynamic range.

# Distributions for count data

Explain when to use common statistical models for counts:

Poisson - normal distribution, large sample number (microarray data that has dynamic range limited maximum due to when the probes max out, and therefore uses the Poisson). Why can we use Poisson for microarray when also there are small sample sizes?

NB - non-normal distribution, low number of biological reps. No max for dynamic range. Don't use technical replicates - biological replicates much more useful (cite paper? or later in md)

# DESeq2 workflow

!image

## Normalization

!image

To make accurate comparisons of abundance between samples, the raw counts, or number of reads aligning to each gene, need to be normalized to account for differences in library depth between samples. If making comparisons between genes within samples, gene length should also be taken into account. 

Different types of normalization methods exist, and two of the most common methods include:
 
 - normalization for library size: comparison between samples
 
 !image
 
 - normalization for gene length: comparison within samples
 
 !image
 
 - normalization for RNA composition: comparison between samples (particularly important for differential expression analyses)
 
 !make_image
 
 *Scaling to library size as a form of normalization makes intuitive sense, given it is expected that sequencing a sample to half the depth will give, on average, half the number of reads mapping to each gene. We believe this is appropriate for normalizing between replicate samples of an RNA population. However, library size scaling is too simple for many biological applications. The number of reads expected to map to a gene is not only dependent on the expression level and length of the gene, but also the composition of the RNA population that is being sampled. Thus, if a large number of genes are unique to, or highly expressed in, one experimental condition, the sequencing 'real estate' available for the remaining genes in that sample is decreased. If not adjusted for, this sampling artifact can force the DE analysis to be skewed towards one experimental condition. Current analysis methods have not accounted for this proportionality property of the data explicitly, potentially giving rise to higher false positive rates and lower power to detect true differences.* [[2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)]

*few highly and differentially expressed genes may have strong influence on the total read count, causing the ratio of total read counts not to be a good estimate for the ratio of expected counts*[[1](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)]
 
 
Several common normalization measures exist to account for these differences:

- **CPM (counts per million):** counts scaled by total number of reads. This measure accounts for sequencing depth only.
- **TPM (transcripts per kilobase million):** counts per length of transcript (kb) per million reads mapped. This measure accounts for both sequencing depth and gene length.
- **RPKM/FPKM (reads/fragments per kilobase of exon per million reads/fragments mapped):** similar to TPM. This measure accounts for both sequencing depth and gene length as well; however, it is **not recommended**.
- **Tool-specific metrics for normalization:** DESeq2 uses a median of ratios method, which accounts for sequencing depth and RNA composition[[1]()]. EdgeR uses a trimmed mean of M values (TMM) method that accounts for sequencing depth, RNA composition, and gene length [[2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)]

### RPKM/FPKM (not recommended)
While TPM and RPKM/FPKM normalization methods both account for sequencing depth and gene length, RPKM/FPKM measures are not recommended. **The reason  is that the normalized count values output by the RPKM/FPKM method are not comparable between samples.** 

If you sum the total number of RPKM/FPKM normalized counts for each sample, the total number will be different between samples. Therefore, you cannot compare the normalized counts for each gene equally between samples. 


**RPKM-normalized counts table**

| gene | sample1 | sample2 |
| ----- |:-----:|:-----:|
| MOV10 | 5.5 | 5.5 |
| ABCD | 73.4 | 21.8 |
| ... | ... | ... |
|Total RPKM-normalized counts | 1,000,000 | 1,500,000 |

For example, using the table above, the MOV10 gene has 5.5 RPKM-normalized counts in sample1 out of 1,000,000 total RPKM-normalized counts in the sample, while sample2 also has 5.5 RPKM-normalized counts associated with MOV10 out of 1,500,000 total RPKM-normalized counts for sample2. Therefore, we cannot directly compare the counts for MOV10 (or any other gene) between sample1 and sample2 because the total number of normalized counts are different between samples. Sample1 has a greater proportion of counts associated with MOV10 (5.5/1,000,000) than does sample2 (5.5/1,500,000) even though the RPKM count values are the same.

### TPM (recommended)
In contrast to RPKM/FPKM, TPM-normalized counts normalize for both sequencing depth and gene length, but have the same total TPM-normalized counts per sample. Therefore, the normalized count values are comparable both between and within samples.

### DESeq2-normalized counts - Median of ratios method
Since tools for differential expression analysis are comparing the counts between sample groups for the same gene, gene length does not need to be accounted for by the tool. However, sequencing depth and RNA composition do need to be included in the normalization.

To normalize for sequencing depth and RNA composition, DESeq2 uses the median of ratios method, which performs the following steps when you run the tool:

**Step 1: creates a pseudo-reference sample (row-wise geometric mean)**
For each gene, a pseudo-reference sample is created that is equal to the geometric mean across all samples.

| gene | sample1 | sample2 | | pseudo-reference sample  |
| ----- |:-----:|:-----:| |:-----:|
| MOV10 | 1489 | 906 | | sqrt(1489 * 906) = 1161.5 |
| ABCD | 24 | 13 | | sqrt(24 * 13) = sqrt(24 * 13) = 17.7 |
| ... | ... | ... | | ... |

**Step 2: calculates ratio of each sample to the reference**

For every gene in a sample, the ratios (sample1/ref) are calculated (as shown below). This is performed for each sample in the dataset.

| gene | sample1 | sample2 | | pseudo-reference sample  | ratio sample1/ref | ratio sample2/ref |
| ----- |:-----:|:-----:| |:-----:| :-----: | :-----: |
| EF2A | 1489 | 906 | | sqrt(1489 * 906) = 1161.5 | 1489/1161.5 = 1.28 | 906/1161.5 = 0.78 |
| BBC1 | 22 | 13 | | sqrt(22 * 13) = 16.9 | 22/16.9 = 1.30| 13/16.9 = 0.77 |
| ... | ... | ... | | ... |

**Step 3: takes sample's median value as that sample's normalization factor**

 The median value of all ratios is taken as the normalization factor (size factor) for that sample. 

<img src="../img/deseq_median_of_ratios.png" width="400">

The median of ratios method makes the assumption that not ALL genes are differentially expressed, and the normalization factors should account for sequencing depth (since the ratios for most genes should be most influenced by sequencing depth) and RNA composition of the sample (large outlier genes will not represent the median ratio values). This method is robust to imbalance in up-/down- regulation and large numbers of differentially expressed genes.



