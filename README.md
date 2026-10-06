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


## 02. post_DESeq2_rhythmicity_Bio_Int.R

This script does some early analysis focusing on 'biological interpretation' of the genes identified as being rhythmic using `DESeq2_Rhythmic_Gene_expression.R`. Specifically, this analysis focuses on the following:

```bash
- Rhythmic gene summary
- Phase distribution
- Phase bins
- Amplitude distribution
- Phase vs amplitude
- Phase-ordered VST heatmap
- g:Profiler GO Biological Process enrichment
- g:Profiler KEGG enrichment
- Overall enrichment
- Phase-specific enrichment
- Diagnostics for recognised/unrecognised gene IDs
```

### Input file
```bash
DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv
```

### Package prerequisites
```R
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(gprofiler2)
```

### Model organism
`g:Profiler` explicitly requires a model organism. Here, we are set to *H. glaber* (`hgfemale`).

### Output files

| Output category                     | File(s)                                        | Explanation                                                                                         |
| ----------------------------------- | ---------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| **Rhythmic gene summary**           | `rhythmic_gene_summary.csv`                    | Overall number and characteristics of rhythmic genes, including amplitude and peak-phase statistics |
| **Genes by phase**                  | `rhythmic_genes_by_phase.csv`                  | Number of rhythmic genes peaking in each 4-hour ZT window                                           |
| **Phase-ordered gene list**         | `rhythmic_genes_final_phase_ordered.csv`       | Complete rhythmic-gene list ordered by peak time                                                    |
| **Top genes per phase**             | `top_10_rhythmic_genes_per_phase.csv`          | Top 10 rhythmic genes within each phase group                                                       |
| **Peak-phase distribution**         | `rhythmic_gene_phase_distribution.pdf/png`     | Visual distribution of rhythmic-gene peak times across the 24-h cycle                               |
| **Phase-bin distribution**          | `rhythmic_genes_by_phase_bin.pdf/png`          | Visual comparison of rhythmic-gene numbers across the six ZT phase bins                             |
| **Amplitude distribution**          | `rhythmic_gene_amplitude_distribution.pdf/png` | Distribution of rhythmicity amplitudes                                                              |
| **Phase vs amplitude**              | `rhythmic_phase_vs_amplitude.pdf/png`          | Relationship between peak timing and rhythm amplitude                                               |
| **VST expression matrix**           | `phase_ordered_rhythmic_gene_VST_matrix.csv`   | VST expression values for rhythmic genes, ordered by peak phase                                     |
| **Expression heatmap**              | `all_rhythmic_genes_phase_ordered_heatmap.pdf` | Visual confirmation of temporal expression patterns across rhythmic genes                           |
| **Overall GO enrichment**           | `GO_Biological_Process_enrichment_all.csv`     | All GO Biological Process terms detected among rhythmic genes                                       |
| **Significant GO enrichment**       | `GO_Biological_Process_enrichment_FDR05.csv`   | GO Biological Processes significantly enriched at **FDR < 0.05**                                    |
| **Overall KEGG enrichment**         | `KEGG_enrichment_all.csv`                      | All KEGG pathways detected among rhythmic genes                                                     |
| **Significant KEGG enrichment**     | `KEGG_enrichment_FDR05.csv`                    | KEGG pathways significantly enriched at **FDR < 0.05**                                              |
| **Overall enrichment plot**         | `gProfiler_overall_enrichment.pdf/png`         | Visual summary of the strongest overall functional-enrichment results                               |
| **Phase-specific GO enrichment**    | `gProfiler_ZT*_GO_BP_all.csv`                  | GO Biological Processes associated with genes peaking in each phase                                 |
| **Phase-specific significant GO**   | `gProfiler_ZT*_GO_BP_FDR05.csv`                | Significant GO processes for each phase (**FDR < 0.05**)                                            |
| **Phase-specific KEGG enrichment**  | `gProfiler_ZT*_KEGG_all.csv`                   | KEGG pathways associated with genes peaking in each phase                                           |
| **Phase-specific significant KEGG** | `gProfiler_ZT*_KEGG_FDR05.csv`                 | Significant KEGG pathways for each phase (**FDR < 0.05**)                                           |
| **Phase enrichment summary**        | `gProfiler_phase_enrichment_combined.csv`      | Combined enrichment results across all six circadian phases                                         |
| **Phase enrichment plots**          | `gProfiler_ZT*_enrichment.pdf/png`             | Visual summaries of significant functional enrichment for each phase                                |

