---
title: "Functional Analysis for RNA-Seq"
author: "Mary Piper"
date: "October 17, 2016"
---

Approximate time: 60 minutes

Learning Objectives:
-------------------

*  Determine how biological functions are attributed to genes using Gene Ontology terms
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
There are a plethora of functional enrichment tools that perform some type of over-representation analysis by querying databases containing information about gene function and interactions. Querying these databases for gene function often requires the use of a **consistent vocabulary** to describe gene function. One of the most widely-used vocabularies is the Gene Ontology (GO). This vocabulary was established by the Gene Ontology project, and the words in the vocabulary are referred to as GO terms. 

### Gene Ontology project

"The Gene Ontology project is a collaborative effort to address the need for consistent descriptions of gene products across databases" [[2](geneontology.org/page/documentation)]. The [Gene Ontology Consortium](http://geneontology.org/page/go-consortium-contributors-list) maintains the GO terms, and these GO terms are incorporated into gene annotations in many of the popular repositories for animal, plant, and microbial genomes. 

Tools that investigate **enrichment of biological functions or interactions** can query these databases for GO terms associated with a list of genes to determine whether any GO terms associated with particular functions or interactions are enriched in the gene set. Therefore, to best use and interpret the results from these functional analysis tools, it is helpful to have a good understanding of the GO terms themselves.

### GO terms

#### GO Ontologies

To describe the roles of genes and gene products, GO terms are organized into three independent controlled vocabularies (ontologies) in a species-independent manner: 

- **Biological process:** refers to the biological process the gene or gene product is involved in (e.g. "transcription", "signal transduction", "apoptosis", etc.). 
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

To determine whether a GO term or pathway is *significantly* over- or under-represented, tools often perform **hypergeometric testing**. "The hypergeometric distribution is a discrete probability distribution that describes the probability of _k_ successes in _n_ draws, without replacement, from a finite population of size _N_ that contains exactly _K_ successes, wherein each draw is either a success or a failure" [[4](https://en.wikipedia.org/wiki/Hypergeometric_distribution)]. 

Therefore, using our example, the hypergeometric distribution describes the probability of 50 genes (k) being associated with "DNA repair", for all genes in our gene list (n=1,000), from a population of all of the genes in entire genome (N=13,000) which contains 85 genes (K) associated with "DNA repair".

The calculation of probability of k successes follows the formula:

![hypergeo](../img/hypergeo.png) 


### gProfiler

[gProfileR](http://biit.cs.ut.ee/gprofiler/index.cgi) is a tool for the interpretation of large gene lists which can be run using a web interface or through R. The core tool takes a gene list as input and performs statistical enrichment analysis using hypergeometric testing to provide interpretation to user-provided gene lists. Multiple sources of functional evidence are considered, including Gene Ontology terms, biological pathways, regulatory motifs of transcription factors and microRNAs, human disease annotations and protein-protein interactions. The user selects the organism and the sources of evidence to test. There are also additional parameters to change various thresholds and tweak the stringency to the desired level. 

![gprofiler](../img/gProfiler.png)

You can use gProfiler for a wide selection of organisms, and the tool accepts your gene list as input. If your gene list is ordered (e.g. by padj. values), then gProfiler will take the order of the genes into account when outputting enriched terms or pathways.

In addition, a large number (70%) of the functional annotations of GO terms are determined using _in silico_ methods to infer function from electronic annotation (IEA). While these annotations can offer valuable information, the information is of lower confidence than experimental and computational studies, and these functional annotations can be easily filtered out. 

The color codes in the gProfiler output represent the quality of the evidence for the functional annotation. For example, weaker evidence is depicted in blue, while strong evidence generated by direct experiment is shown with red or orange. Similar coloring is used for pathway information, with well-researched pathway information shown in black, opposed to lighter colors. Grey coloring suggests an unknown gene product or annotation. For more information, please see the [gProfiler paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933153/).

Also, due to the hierarchical structure of GO terms, you may return many terms that seem redundant since they are child and parent terms. gProfiler allows for 'hierarchical filtering', returning only the best term per parent term.

#### Running gProfiler

For our gProfiler analysis, we are going to use the Mov10_oe results table subsetted using the thresholds of logFC = 1 (2x) and padj = 0.05, which is saved as `SigOE`. 
```
### Functional analysis of MOV10 Overexpression using gProfileR (some of these are defaults; check help pages) 

library(gProfileR)

gprofiler_results_oe <- gprofiler(query = sigOE, 
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

write.table(gprofiler_results_oe, 
            "results/gprofiler_MOV10_oe.txt", 
            sep="\t", quote=F, row.names=F)
```

Now, extract only the lines in the gProfiler results with GO term accession numbers for downstream analyses:

```
## Extract GO IDs for downstream analysis

allterms_oe <- gprofiler_results_oe$term.id

GOs_oe <- allterms_oe[grep('GO:', allterms_oe)]

write.table(GOs_oe, "results/GOs_oe.txt", sep="\t", quote=F, row.names=F, col.names=F)
```

### REVIGO

[REVIGO](http://revigo.irb.hr/) is a web-based tool that can take our list of GO terms, collapse redundant terms by semantic similarity, and summarize them graphically. 

![REVIGO_input](../img/revigo_input.png)

Open `GOs_oe.txt` and copy and paste the GO ids into the REVIGO search box, and submit.

![REVIGO_output](../img/revigo_output.png)

***gProfiler and REVIGO are great tools to validate experimental results and to make hypotheses. These tools suggest pathways that may be involved with your condition of interest; you should NOT use these tools to make conclusions about the pathways involved in your experimental process.***

J. Reimand, T. Arak, P. Adler, L. Kolberg, S. Reisberg, H. Peterson, J. Vilo. g:Profiler -- a web server for functional interpretation of gene lists (2016 update). Nucleic Acids Research 2016; doi: 10.1093/nar/gkw199

Supek F, Bošnjak M, Škunca N, Šmuc T. REVIGO summarizes and visualizes long lists of Gene Ontology terms. PLoS ONE 2011. doi:10.1371/journal.pone.0021800

## [Other functional analysis methods](https://github.com/hbc/DGE_workshop/blob/master/lessons/functional_analysis_other_methods.md)

Over-representation analyses are only a single type of functional analysis method that is available for teasing apart the biological processes important to your condition of interest. Other types of analyses can be equally important or informative, including functional class scoring and pathway topology methods. Functional class scoring methods most often take as input the foldchanges for all genes, then look to see whether gene sets for particular biological processes are enriched among the high or low fold changes. This type of analysis can be particularly helpful if differential expression analysis only output a small list of significant DE genes. Finally, pathway topology analysis often takes into account both fold changes and adjusted p-values to identify dysregulated pathways and outputs whether pathways are inhibited/activated. We have [materials](https://github.com/hbc/DGE_workshop/blob/master/lessons/functional_analysis_other_methods.md) to lead you through these other types of functional analyses, and we encourage you to take the time to work through them.

![Pathway analysis tools](../img/pathway_analysis.png)

## Resources for functional analysis

* g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
* DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
* clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
* GeneMANIA - http://www.genemania.org/
* GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
* WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
* AmiGO - http://amigo.geneontology.org/amigo
* ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
* WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
* GSEA - http://software.broadinstitute.org/gsea/index.jsp
* SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
* Gage/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
