# cosinor_expression_analysis
A series of (predominantly) R scripts to explore circadian rhythms in time-series transcriptomic sequence data.

## 01. DESeq2_Rhythmic_Gene_expression.R

A cosinor analytical approach to differential expression analysis using [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 

Assumes the following (which are adjustable):

* **Experiment**
```bash
6 timepoints (ZT0, ZT4, ZT8, ZT12, ZT16, ZT20)
4 biological replicates per time point
```
* **Model**
```bash
Full    = ~ cos24 + sin24
Reduced = ~ 1
```
* **Statistical test**
```bash
DESeq2 likelihood-ratio test (LRT)
```
* **Cosinor plots**
```bash
ONLY the TOP_N genes with FDR < 0.05
Ranked by lowest adjusted p-value
```

### Package prerequisites
```R
library(DESeq2)
library(ggplot2)
library(pheatmap)
```

### Output files
| Output                                               | Explanation                                                                     |
| ---------------------------------------------------- | ------------------------------------------------------------------------------- |
| `DESeq2_24h_rhythm_all_genes.csv`                    | DESeq2 LRT results for all genes tested                                         |
| `DESeq2_24h_rhythmic_genes_FDR05.csv`                | Genes with **FDR/padj < 0.05** — your significant rhythmic genes                |
| `DESeq2_24h_rhythm_with_phase_amplitude.csv`         | All genes plus cosine coefficient, sine coefficient, amplitude, and peak ZT     |
| `DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv` | Same information restricted to **FDR < 0.05** genes                             |
| `top_50_rhythmic_genes_by_FDR.csv`                   | The 50 most significant rhythmic genes, ranked by lowest FDR                    |
| `DESeq2_normalized_counts.csv`                       | DESeq2 size-factor-normalized expression values                                 |
| `DESeq2_VST_expression.csv`                          | Variance-stabilized expression values                                           |
| `PCA_ZT.pdf`                                         | PCA showing sample-level expression structure                                   |
| `sample_correlation_heatmap.pdf`                     | Correlation between biological replicates/samples                               |
| `DESeq2_rhythm_MAplot.pdf`                           | MA plot of the rhythmicity analysis                                             |
| `top_20_rhythmic_genes_heatmap.pdf`                  | Heatmap of the top 20 FDR-significant rhythmic genes                            |
| `top_50_rhythmic_genes_cosinor_plots.pdf`            | Combined PDF containing cosinor plots for only the top 50 genes                 |
| `top_50_cosinor_plots/`                              | Individual PDF and PNG cosinor plots for those 50 genes                         |
| `DESeq2_24h_rhythm_dds.rds`                          | Complete DESeq2 object for later analysis                                       |


## 02. 
