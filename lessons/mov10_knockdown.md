## Define contrasts, extract results table and shrink log2 fold changes

```{r}
contrast <-  c("sampletype", "MOV10_knockdown", "control")

res_tableKD <- results(dds, contrast=contrast)

res_tableKD <- lfcShrink(dds, contrast=contrast, res=res_tableKD)
```

## Set a threshold and subset the results

```{r}
threshold_KD <- res_tableKD$padj < padj.cutoff & abs(res_tableKD$log2FoldChange) > lfc.cutoff

res_tableKD$threshold <- threshold_KD

sigKD <- subset(res_tableKD, threshold)
```
