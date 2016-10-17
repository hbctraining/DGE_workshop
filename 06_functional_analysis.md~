---
title: "Functional Analysis for RNA-Seq"
author: "Mary Piper"
date: "Tuesday, July 5th, 2016"
---

Approximate time: 105 minutes

Learning Objectives:
-------------------

*  Determine how functions are attributed to genes using Gene Ontology terms
*  Understand the theory of how functional enrichment tools yield statistically enriched functions or interactions
*  Discuss functional analysis using over-representation analysis, functional class scoring, and pathway topology methods
*  Explore functional analysis tools

# Functional analysis 

The output of RNA-Seq differential expression analysis is a list of significant differentially expressed genes (DEGs). To gain greater biological insight on the DEGs there are various analyses that can be done:

- determine whether there is enrichment of known biological functions, interactions, or pathways
- identify genes' involvement in novel pathways or networks by grouping genes together based on similar trends
- use global changes in gene expression by visualizing all genes being significantly up- or down-regulated in the context of external interaction data

Generally for any differential expression analysis, it is useful to interpret the resulting gene lists using freely available web- and R-based tools.  While tools for functional analysis span a wide variety of techniques, they can loosely be categorized into three main types: over-representation analysis, functional class scoring, and pathway topology [[1](../../resources/pathway_tools.pdf)]. 

![Pathway analysis tools](../img/pathway_analysis.png)

## Over-representation analysis
There are a plethora of functional enrichment tools that perform some type of over-representation analysis by querying databases containing information about gene function and interactions. Querying these databases for gene function requires the use of a **consistent vocabulary** to describe gene function. One of the most widely-used vocabularies is the Gene Ontology (GO). This vocabulary was established by the Gene Ontology project, and the words in the vocabulary are referred to as GO terms. 

### Gene Ontology project

"The Gene Ontology project is a collaborative effort to address the need for consistent descriptions of gene products across databases" [[2](geneontology.org/page/documentation)]. The [Gene Ontology Consortium](http://geneontology.org/page/go-consortium-contributors-list) maintains the GO terms, and these GO terms are incorporated into gene annotations in many of the popular repositories for animal, plant, and microbial genomes. 

Tools that investigate **enrichment of biological functions or interactions** can query these databases for GO terms associated with a list of genes to determine whether any GO terms associated with particular functions or interactions are enriched in the gene set. Therefore, to best use and interpret the results from these functional analysis tools, it is helpful to have a good understanding of the GO terms themselves.

### GO terms

#### GO Ontologies

To describe the roles of genes and gene products, GO terms are organized into three independent controlled vocabularies (ontologies) in a species-independent manner: 

- **Biological process:** refers to the biological role involving the gene or gene product, and could include "transcription", "signal transduction", and "apoptosis". A biological process generally involves a chemical or physical change of the starting material or input.
- **Molecular function:** represents the biochemical activity of the gene product, such activities could include "ligand", "GTPase", and "transporter". 
- **Cellular component:** refers to the location in the cell of the gene product. Cellular components could include "nucleus", "lysosome", and "plasma membrane".

"The relationships between a gene product to biological process, molecular function and cellular component are one-to-many, reflecting the biological reality that a particular protein may function in several processes, contain domains that carry out diverse molecular functions, and participate in multiple alternative interactions with other proteins, organelles or locations in the cell" [[3](go.pdf)]. Therefore, a single gene product can be associated with many GO terms. Each GO term has a term name (e.g. DNA repair) and a unique term accession number (GO:0005125).

#### GO term hierarchy

Some gene products are well-researched, with vast quantities of data available regarding their biological processes and functions. However, other gene products have very little data available about their roles in the cell. 

For example, the protein, "p53", would contain a wealth of information on it's roles in the cell, whereas another protein might only be known as a "membrane-bound protein" with no other information available. 

The GO ontologies were developed to describe and query biological knowledge with differing levels of information available. To do this, GO ontologies are loosely hierarchical, ranging from general, 'parent', terms to more specific, 'child' terms. The GO ontologies are "loosely" hierarchical since 'child' terms can have multiple 'parent' terms.

Some genes with less information may only be associated with general 'parent' terms or no terms at all, while other genes with a lot of information have many terms.

![Nature Reviews Cancer 7, 23-34 (January 2007)](../img/go_heirarchy.jpg)

[Tips for working with GO terms](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003343)

### Hypergeometric testing

In a set of genes, the frequency of GO terms can be determined, and the comparison of frequencies between a gene list & a “background” set will inform us about the over- or under- representation of the GO terms. This type of testing can inform us about over- or under-representation of other entities such as particular motifs or pathways too.

![go_frequencies](../img/go_freq.png)

To determine whether GO terms (or motifs and pathways) are over- or under-represented, you can determine the probability of having a certain number of genes associated with specific GO terms for the size of the gene list based on the background set. The background dataset can be all genes in genome for your organism or you can select your own background to use.

For example, let's suppose there are 13,000 total genes in the honeybee genome and 85 genes are associated with the GO term "DNA repair". In your gene list, there are 50 genes associated with "DNA repair" out of 1,000 genes in gene list. 

By comparing the ratios, 85/13,000 in "background" dataset and 50/1,000 in your gene list, it's evident that the GO term "DNA repair" is over-represented in your dataset.

To determine whether a GO term or pathway is significantly over- or under-represented, tools often perform **hypergeometric testing**. "The hypergeometric distribution is a discrete probability distribution that describes the probability of _k_ successes in _n_ draws, without replacement, from a finite population of size _N_ that contains exactly _K_ successes, wherein each draw is either a success or a failure" [[4](https://en.wikipedia.org/wiki/Hypergeometric_distribution)]. 

Therefore, using our example, the hypergeometric distribution describes the probability of 50 genes (k) being associated with "DNA repair", for all genes in our gene list (n=1,000), from a population of all of the genes in entire genome (N=13,000) which contains 85 genes (K) associated with "DNA repair".

The calculation of probability of k successes follows the formula:

![hypergeo](../img/hypergeo.png) 


### gProfiler

[gProfileR](http://biit.cs.ut.ee/gprofiler/index.cgi) is a web-based tool for the interpretation of large gene lists. The core tool takes a gene list as input and performs statistical enrichment analysis using hypergeometric testing to provide interpretation to user-provided gene lists. Multiple sources of functional evidence are considered, including Gene Ontology terms, biological pathways, regulatory motifs of transcription factors and microRNAs, human disease annotations and protein-protein interactions. The user selects the organism and the sources of evidence to test. There are also additional parameters to change various thresholds and tweak the stringency to the desired level. 

![gprofiler](../img/gProfiler.png)

You can use gProfiler for a wide selection of organisms, and the tool accepts your gene list as input. If your gene list is ordered (e.g. by padj. values), then gProfiler will take the order of the genes into account when outputting enriched terms or pathways.

In addition, a large number (70%) of the functional annotations of GO terms are determined using _in silico_ methods to infer function from electronic annotation (IEA). While these annotations can offer valuable information, the information is of lower confidence than experimental and computational studies, and these functional annotations can be easily filtered out. 

The color codes in the gProfiler output represent the quality of the evidence for the functional annotation. For example, weaker evidence is depicted in blue, while strong evidence generated by direct experiment is shown with red or orange. Similar coloring is used for pathway information, with well-researched pathway information shown in black, opposed to lighter colors. Grey coloring suggests an unknown gene product or annotation. For more information, please see the [gProfiler_paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933153/).

Also, due to the hierarchical structure of GO terms, you may return many terms that seem redundant since they are child and parent terms. gProfiler allows for 'hierarchical filtering', returning only the best term per parent term.

#### Running gProfiler

For our gProfiler analysis, we are going to use the Mov10_kd results table subsetted using the thresholds of logFC = 0.58 (1.5x) and padj = 0.05. 
**You can download the subsetted Mov10_kd results file using this [link](https://github.com/hbc/NGS_Data_Analysis_Summer2016/raw/master/sessionIII/data/Mov10_kd_logFC_0.58_pVal_0.05.txt).**

Save the file to your results directory in your `DEanalysis` project.


* Under **Options**: keep all defaults checked and for `Hierarchical Filtering` use the pulldown to select `Best per parent`
* Click on **Advanced Options** and for `Significance Threshold` select `Benjamini-Hochberg FDR`
* From the functional evidence selections choose the following: Gene Ontology (biological process, molecular function), [KEGG](http://www.genome.jp/kegg/) ([KEGG paper](http://nar.oxfordjournals.org/content/44/D1/D457.full.pdf)), and [Reactome](http://www.reactome.org).
* Press **g:Profile!** 


> Take a look at the list of terms that appear. Do you see anything relevant, given what you know about this dataset? Run the analysis again but this time change the appropriate parameter to export your results to file. 

#### gProfiler in R

While the web interface for gProfiler is a bit more intuitive to understand, we don't actually need to leave R to run gProfiler. 

Go to RStudio and click on the File menu and select 'Open project ...'

Navigate to `~/Desktop/DEanalysis/` and double click on the `DEanalysis.Rproj` file.

```
### Functional analysis of MOV10 Overexpression using gProfileR (some of these are defaults; check help pages) 

library(gProfileR)

sigKD_FC58 <- scan(file="results/Mov10_kd_logFC_0.58_pVal_0.05.txt", what="character")

gprofiler_results_kd <- gprofiler(query = sigKD_FC58, 
                                  organism = "hsapiens",
                                  ordered_query = F, 
                                  exclude_iea = F, 
                                  max_p_value = 0.05, 
                                  max_set_size = 0,
                                  correction_method = "fdr",
                                  hier_filtering = "none", 
                                  domain_size = "annotated",
                                  custom_bg = "")

```

Let's save the gProfiler results to file:

```
## Write results to file

write.table(gprofiler_results_kd, 
            "results/gprofiler_MOV10_kd.txt", 
            sep="\t", quote=F, row.names=F)
```

Now, extract only the lines in the gProfiler results with GO term accession numbers for downstream analyses:

```
## Extract GO IDs for downstream analysis

allterms_kd <- gprofiler_results_kd$term.id

GOs_kd <- allterms_kd[grep('GO:', allterms_kd)]

write.table(GOs_kd, "results/GOs_kd.txt", sep="\t", quote=F, row.names=F, col.names=F)
```

### REVIGO

[REVIGO](http://revigo.irb.hr/) is a web-based tool that can take our list of GO terms, collapse redundant terms by semantic similarity, and summarize them graphically. 


![REVIGO_input](../img/revigo_input.png)

Open `GOs_oe.txt` and copy and paste the GO ids into the REVIGO search box, and submit.

![REVIGO_output](../img/revigo_output.png)

***gProfiler and REVIGO are great tools to validate experimental results and to make hypotheses. These tools suggest pathways that may be involved with your condition of interest, and you should NOT use these tools to make conclusions about the pathways involved in your experimental process.***

J. Reimand, T. Arak, P. Adler, L. Kolberg, S. Reisberg, H. Peterson, J. Vilo. g:Profiler -- a web server for functional interpretation of gene lists (2016 update). Nucleic Acids Research 2016; doi: 10.1093/nar/gkw199

Supek F, Bošnjak M, Škunca N, Šmuc T. REVIGO summarizes and visualizes long lists of Gene Ontology terms. PLoS ONE 2011. doi:10.1371/journal.pone.0021800

## Functional class scoring tools
Functional class scoring (FCS) tools, such as [GSEA](http://software.broadinstitute.org/gsea/index.jsp), use the gene-level statistics from the differential expression results to determine pathway-level expression changes. The hypothesis of FCS methods is that although large changes in individual genes can have significant effects on pathways (and will be detected via ORA methods), weaker but coordinated changes in sets of functionally related genes (i.e., pathways) can also have significant effects.  Thus, rather than setting an arbitrary threshold to identify 'significant genes', **all genes are considered** in the analysis. The gene-level statistics from the dataset are aggregated to generate a single pathway-level statistic and statistical significance of each pathway is reported.

### Gene set enrichment analysis using GAGE and Pathview
Using the log2 fold changes obtained from the DESeq2 analysis for every gene, gene set enrichment analysis and pathway analysis was performed using [GAGE (Generally Applicable Gene-set Enrichment for Pathway Analysis)](http://bioconductor.org/packages/release/bioc/html/gage.html) and [Pathview](http://bioconductor.org/packages/release/bioc/html/pathview.html) tools.

For gene set or pathway analysis using GAGE, coordinated differential expression over gene sets is tested instead of changes of individual genes. "Gene sets are pre-defined groups of genes, which are functionally related. Commonly used gene sets include those derived from KEGG pathways, Gene Ontology terms, gene groups that share some other functional annotations, etc. Consistent perturbations over such gene sets frequently suggest mechanistic changes." [[1](https://www.bioconductor.org/packages/devel/bioc/vignettes/gage/inst/doc/gage.pdf)]

"GAGE assumes a gene set comes from a different distribution than the background and uses two-sample t-test to account for the gene set specific variance as well as the background variance. The two-sample t-test used by GAGE identifies gene sets with modest but consistent changes in gene expression level."[[2](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161)]

Pathview allows for the integration of the data generated by GAGE and visualization of the pathways from the dataset.


#### Exploring enrichment of KEGG pathways
To get started with GAGE and Pathview analysis, we need to load multiple libraries:

```
# Loading the packages needed for GAGE and Pathview analysis

library(gage)
library(pathview)
library(gageData)
library(biomaRt)
library(org.Hs.eg.db)
```

To determine whether pathways in our dataset are enriched, we need to first obtain the gene sets to test:

```
# Create datasets with KEGG gene sets to test

kegg_human <- kegg.gsets(species = "human", id.type = "kegg")
names(kegg_human)

kegg.gs <- kegg_human$kg.sets[kegg_human$sigmet.idx]
head(kegg.gs)
```

Now that we have our pathways to test, we need to bring in our own data. We will use the log2 fold changes output by DESeq2, `res_tableKD_sorted`, to determine whether particular pathways are enriched. GAGE requires the genes have Entrez IDs, so we will also convert our gene names into Entrez IDs prior to analysis. *If you do not have this file, you can download it using this [link](https://github.com/hbc/NGS_Data_Analysis_Summer2016/raw/master/sessionIII/data/res_tableKD.txt).*

> A useful tutorial for using GAGE and Pathview is available from Stephen Turner on R-bloggers: [http://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/](http://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/)

```
#Set up

## Turn DESeq2 results into a dataframe

DEG <- data.frame(res_tableKD_sorted)

## Tool expects entrez IDs, so use Biomart to acquire IDs for the gene names

### Create database to search

mart <- useDataset("hsapiens_gene_ensembl", 
                   useMart('ENSEMBL_MART_ENSEMBL', 
                           host =  'www.ensembl.org')) 
                           
###List entrez IDs for each gene name

entrez <- getBM(filters= "external_gene_name", 
                attributes= c("external_gene_name", "entrezgene"),
                values= row.names(DEG),
                mart= mart)
## Create a column in results dataset, DEG, to merge with the DEG dataset with entrez IDs

DEG$external_gene_name <- row.names(DEG)

## Merge results dataset, DEG, with Entrez IDs dataset

entrez_results <- merge(DEG, entrez, by="external_gene_name")
head(entrez_results, n=15)

entrez_results <- subset(entrez_results, entrezgene != "NA")
head(entrez_results, n=15)

## Extract only the log2FC values

foldchanges <- entrez_results$log2FoldChange
head(foldchanges)

## Name foldchanges vector with corresponding Entrez IDs as needed by GAGE tool

names(foldchanges) <- entrez_results$entrezgene
head(foldchanges)
```
Now the data is ready to run for GAGE analysis. The [GAGE vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/gage/inst/doc/gage.pdf) provides detailed information on running the analysis. Note that you can run the analysis and look for pathways with genes statistically only up- or down-regulated. Alternatively, you can explore statistically perturbed pathways, which are enriched in genes that may be either up- or down- regulated. For KEGG pathways, looking at both types of pathways could be useful.

```
# Run GAGE

keggres = gage(foldchanges, gsets=kegg.gs, same.dir=T)
names(keggres)
head(keggres$greater) #Pathways that are up-regulated
head(keggres$less) #Pathways that are down-regulated

# Explore genes that are up-regulated

sel_up <- keggres$greater[, "q.val"] < 0.05 & !is.na(keggres$greater[, "q.val"])
path_ids_up <- rownames(keggres$greater)[sel_up]
path_ids_up

# Get the pathway IDs for the significantly up-regulated pathways

keggresids = substr(path_ids_up, start=1, stop=8)
keggresids
```
Now that we have the IDs for the pathways that are significantly up-regulated in our dataset, we can visualize these pathways and the genes identified from our dataset causing these pathways to be enriched using [Pathview](https://www.bioconductor.org/packages/devel/bioc/vignettes/pathview/inst/doc/pathview.pdf). 

```
# Run Pathview

## Use Pathview to view significant up-regulated pathways

pathview(gene.data = foldchanges, pathway.id=keggresids, species="human", kegg.dir="results/")
```
![spliceosome](../img/hsa03040.pathview.png)

#### Exploring enrichment of biological processes using GO terms in GAGE
Using the GAGE tool, we can identify significantly enriched gene ontology terms for biological process and molecular function based on the log2 fold changes for all genes. While gProfileR is an overlap statistic analysis tool which uses a threshold (adjusted p<0.05 here) to define which genes are analyzed for GO enrichment, gene set enrichment analysis tools like GAGE use a list of genes (here ranked by logFC) without using a threshold. This allows GAGE to use more information to identify enriched biological processes. The introduction to GSEA goes into more detail about the advantages of this approach: [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/).

```
#Acquire datasets

data(go.sets.hs)
head(names(go.sets.hs))

data(go.subs.hs)
names(go.subs.hs)
head(go.subs.hs$MF)

#Use gage to explore enriched biological processes
#Biological process 

go_bp_sets = go.sets.hs[go.subs.hs$BP]
```

> If we wanted to identify enriched molecular functions we would use the code: `go.sets.hs[go.subs.hs$MF]`


```
# Run GAGE
go_bp_res = gage(foldchanges, gsets=go_bp_sets, same.dir=T)
class(go_bp_res)
names(go_bp_res)
head(go_bp_res$greater)
go_df_enriched <- data.frame(go_bp_res$greater)

GO_enriched_BP <- subset(go_df_enriched, q.val < 0.05)
GO_enriched_BP

write.table(GO_enriched_BP, "Mov10_GAGE_GO_BP.txt", quote=F)
```

Weijun Luo, Michael Friedman, Kerby Shedden, Kurt Hankenson, and Peter Woolf. GAGE: generally applicable
gene set enrichment for pathway analysis. BMC Bioinformatics, 2009. doi:10.1186/1471-2105-10-161.

Weijun Luo and Cory Brouwer. Pathview: an R/Bioconductor package for pathway-based data integration
and visualization. Bioinformatics, 29(14):1830-1831, 2013. doi: 10.1093/bioinformatics/btt285.


## Pathway topology tools
The previous analyses did not explore how genes interact with each other (e.g. activation, inhibition, phosphorylation, ubiquitination, etc) to determine the pathway-level statistics. Pathway topology-based methods utilize the number and type of interactions between gene product (our DE genes) and other gene products to infer gene function or pathway association. 

### SPIA
The [SPIA (Signaling Pathway Impact Analysis)](http://bioconductor.org/packages/release/bioc/html/SPIA.html) tool can be used to integrate the lists of differentially expressed genes determined by DESeq2, their fold changes, and pathway topology to identify affected pathways. The blog post from [Getting Genetics Done](http://www.gettinggeneticsdone.com/2012/03/pathway-analysis-for-high-throughput.html) provides a step-by-step procedure for using and understanding SPIA.


Before we run SPIA, we need to remove all NA values and duplicated Entrez IDs:

```
# Set-up

## Significant genes is a vector of fold changes where the names are ENTREZ gene IDs. The background set is a vector of all the genes represented on the platform.

## Convert ensembl to entrez ids

entrez_results <- merge(DEG, entrez, by="external_gene_name")

sig_genes <- subset(entrez_results, padj< 0.05)$log2FoldChange

names(sig_genes) <- subset(entrez_results, padj< 0.05)$entrezgene

head(sig_genes)


## Remove NA and duplicated values
sig_genes <- sig_genes[!is.na(names(sig_genes))] 

sig_genes <- sig_genes[!duplicated(names(sig_genes))]

background_genes <- entrez_results$entrezgene

background_genes <- background_genes[!duplicated(background_genes)]

```

Now that we have our background and significant genes in the appropriate format, we can run SPIA:

```
# Run SPIA.

library(SPIA)

spia_result <- spia(de=sig_genes, all=background_genes, organism="hsa")

head(spia_result, n=20)
```

SPIA outputs a table showing significantly dysregulated pathways based on over-representation and signaling perturbations accumulation. The table shows the following information: `pSize` is the number of genes on the pathway; `NDE` is the number of DE genes per pathway; `tA` is the observed total preturbation accumulation in the pathway; `pNDE` is the probability to observe at least NDE genes on the pathway using a hypergeometric model; `pPERT` is the probability to observe a total accumulation more extreme than tA only by chance; `pG` is the p-value obtained by combining pNDE and pPERT; `pGFdr` and `pGFWER` are the False Discovery Rate and respectively Bonferroni adjusted global p-values; and the Status gives the direction in which the pathway is perturbed (activated or inhibited). KEGGLINK gives a web link to the KEGG website that displays the pathway image with the differentially expressed genes highlighted in red.

We can view the significantly dysregulated pathways by viewing the over-representation and perturbations for each pathway.

```
plotP(spia_result,threshold=0.05)
```
In this plot each pathway is a point and the coordinates are the log of pNDE (using a hypergeometric model) and the p-value from perturbations, pPERT. The oblique lines in the plot show the significance regions based on the combined evidence.

If we choose to explore the significant genes from our dataset occurring in these pathways, we can subset our SPIA results:

```
## Look at pathway 05222 and view kegglink
subset(spia_result, ID == "05222")
```

Then, click on the KEGGLINK, we can view the genes within our dataset from these perturbed pathways:
![perturbed_pathway](../img/hsa05222.png)

Tarca AL, Kathri P and Draghici S (2013). SPIA: Signaling Pathway Impact Analysis (SPIA) using combined evidence of pathway over-representation and unusual signaling perturbations. [http://bioinformatics.oxfordjournals.org/cgi/reprint/btn577v1](http://bioinformatics.oxfordjournals.org/cgi/reprint/btn577v1).

***

## Other Tools

### GeneMANIA

[GeneMANIA](http://genemania.org/) is a tool for predicting the function of your genes. Rather than looking for enrichment, the query gene set is evaluated in the context of curated functional association data and results are displayed in the form of a network. Association data include protein and genetic interactions, pathways, co-expression, co-localization and protein domain similarity. Genes are represented as the nodes of the network and edges are formed by known association evidence. The query gene set is highlighted and so you can find other genes that are related based on the toplogy in the network. This tool is more useful for smaller gene sets (< 400 genes), as you can see in the figure below our input results in a bit of a hairball that is hard to interpret.

![genemania](../img/genemania.png)

> Use the significant gene list generated from the analysis we performed in class as input to GeneMANIA. Using only pathway and coexpression data as evidence, take a look at the network that results. Can you predict anything functionally from this set of genes? 

### Co-expression clustering

Co-expression clustering is often used to identify genes of novel pathways or networks by grouping genes together based on similar trends in expression. These tools are useful in identifying genes in a pathway, when their participation in a pathway and/or the pathway itself is unknown. These tools cluster genes with similar expression patterns to create 'modules' of co-expressed genes which often reflect functionally similar groups of genes. These 'modules' can then be compared across conditions or in a time-course experiment to identify any biologically relevant pathway or network information.

You can visualize co-expression clustering using heatmaps, which should be viewed as suggestive only; serious classification of genes needs better methods.  

The way the tools perform clustering is by taking the entire expression matrix and computing pair-wise co-expression values. A network is then generated from which we explore the topology to make inferences on gene co-regulation. The [WGCNA](http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork ) package (in R) is one example of a more sophisticated method for co-expression clustering.

## Resources for functional analysis

* g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
* DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
* GeneMANIA - http://www.genemania.org/
* GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
* WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
* AmiGO - http://amigo.geneontology.org/amigo
* ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
* WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
* GSEA - http://software.broadinstitute.org/gsea/index.jsp

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
