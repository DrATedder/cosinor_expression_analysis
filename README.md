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


## 03. clock_target_genes.R

Following on from `DESeq2_Rhythmic_Gene_expression.R` (i.e. must be run first), `clock_target_genes.R` pulls out key information about expression and rhythmicity in well known [core mammalian clock genes](https://en.wikipedia.org/wiki/Circadian_clock) (see below).

### Clock genes covered


|Target |Geneid|All_matches     |
|-------|------|----------------|
|CLOCK  |Clock |Clock           |
|BMAL1  |Bmal1 |Bmal1           |
|PER1   |Per1  |Per1            |
|PER2   |Per2  |Per2            |
|PER3   |Per3  |Per3            |
|CRY1   |Cry1  |Cry1            |
|CRY2   |Cry2  |Cry2            |
|REV-ERB|Nr1d2 |Nr1d2; Nr1d1    |
|ROR    |Rorb  |Rorb; Rora; Rorc|
|DBP    |Dbp   |Dbp             |
|NFIL3  |Nfil3 |Nfil3           |

### Package prerequisites

```R
library(ggplot2)
library(patchwork)
```


### Output files

| File(s)                                      | Format | Explanation                                                                                                                                                                                                             |
| -------------------------------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `target_gene_matches.csv`                    | CSV    | Maps the 11 target clock genes to matching gene identifiers in the count matrix. Records the first matching `Geneid` and all matches found using the specified aliases. Unmatched genes are recorded as `NA`.            |
| `target_clock_gene_raw_counts.csv`           | CSV    | Contains the raw, unnormalised read counts for each identified target gene across all samples. Rows are labelled with the target gene names.                                                                             |
| `target_clock_gene_raw_counts_labelled.csv`  | CSV    | Similar to the raw counts table, but includes both the target gene name and its original `Geneid` as explicit columns.                                                                                                   |
| `target_clock_gene_raw_counts_long.csv`      | CSV    | Stores raw counts in long format, with one row per gene–sample combination. Includes the target gene name, original gene identifier, sample name and raw count. Useful for downstream plotting and statistical analysis. |
| `CLOCK_rhythm.pdf`, `BMAL1_rhythm.pdf`, etc. | PDF    | Individual plots showing the normalised expression pattern across Zeitgeber time (ZT) for each identified target gene. Suitable for publication and detailed inspection.                                                 |
| `CLOCK_rhythm.png`, `BMAL1_rhythm.png`, etc. | PNG    | High-resolution (300 dpi) versions of the individual gene rhythm plots, suitable for presentations, reports and image-based documents.                                                                                   |
| `Figure1_Core_Clock_Genes.pdf`               | PDF    | Combines the available plots for six core clock genes: CLOCK, BMAL1, PER1, PER2, PER3 and CRY1. Arranged in a two-column layout.                                                                                         |
| `Figure1_Core_Clock_Genes.png`               | PNG    | High-resolution (300 dpi) version of Figure 1, containing the available core clock gene plots.                                                                                                                           |
| `Figure2_Additional_Clock_Genes.pdf`         | PDF    | Combines the available plots for five additional clock genes: CRY2, REV-ERB, ROR, DBP and NFIL3. Arranged in a two-column layout.                                                                                        |
| `Figure2_Additional_Clock_Genes.png`         | PNG    | High-resolution (300 dpi) version of Figure 2, containing the available additional clock gene plots.                                                                                                                     |
| `target_clock_gene_summary.csv`              | CSV    | Provides a summary of the rhythmic expression statistics for each identified target gene, including the intercept, cosine and sine coefficients, amplitude, peak ZT and FDR-adjusted p-value.                            |
| `target_matches.rds`                         | RDS    | Saves the complete target gene matching results as an R object for reuse in R without repeating the lookup.                                                                                                              |
| `target_plots.rds`                           | RDS    | Saves the individual gene plots as a named R list, allowing them to be reloaded and recombined or modified later.                                                                                                        |
| `target_summary.rds`                         | RDS    | Saves the target gene summary table as an R object for further analysis in R.                                                                                                                                            |

## 0.4. tissue_specific_heatmap.py

This script (**note** you need to actually pull the required genes from your `featureCounts` or equivalent, first), `log₂(count + 1)` normalises the expression data, and produces a simple heatmap to allow visual assessment of expression across samples for key tissue specific gene expression candidates.

### Tissue specific gene expression candidates

```python
hypothalamic_genes = [
    "Agrp", "Npy", "Pomc", "Sst", "Sim1", "Avp", "Oxt",
    "Pdyn", "Ghrh", "Crh", "Pmch", "Kiss1", "Trh",
    "Otp", "Hcrt", "Cartpt"
]

# 
mouse_model_genes = [
    "Fezf1", "Gal", "Gabrq", "Slc18a2", "Magel2",
    "Slc6a3", "Gpx3", "Ngb", "Baiap3"
]

# low or 'negligible' expression in hypothalamus tissue
negative_marker_genes = [
    "Alb", "Cpa1", "Krt1", "Apoa1", "Tnnt2",
    "Krt10", "Slc4a1"
]

# Non-hypothalamic genes expressed in nearby brain tissue
nearby_brain_genes = [
    "Socs6", "Rab37", "Gbx2", "Syt9", "Amotl1",
    "Vangl1", "Prkcd", "Ptpn3", "Tcf7l2", "Slitrk6",
    "Plekhg1", "Rgs16", "Ramp3", "Lef1", "Synpo2",
    "Tnnt1", "Gjc1"
]

```
