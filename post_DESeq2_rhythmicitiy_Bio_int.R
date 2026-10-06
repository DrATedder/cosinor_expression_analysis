############################################################
# DESeq2 24-h rhythmicity: biological interpretation
#
# INPUT:
#   DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv
#
# IMPORTANT:
#   The FIRST COLUMN of the input CSV contains gene IDs.
#
# BACKGROUND:
#   DESeq2_24h_rhythm_all_genes.csv
#
# VST:
#   DESeq2_VST_expression.csv
#
# ORGANISM:
#   NMR = Hglaber
#
# ENRICHMENT:
#   g:Profiler
#
# Databases:
#   GO Biological Process
#   KEGG
#   Reactome
#
# IMPORTANT CHANGE:
#   g:Profiler is run with significant = FALSE.
#
#   This means ALL tested enrichment terms are returned.
#
#   We then separately export:
#
#     *_enrichment_all.csv
#
#     *_enrichment_FDR05.csv
#
#   This prevents an empty KEGG file from hiding the fact
#   that KEGG pathways may have been tested but did not
#   survive FDR correction.
############################################################


############################################################
# 1. PACKAGES
############################################################

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "pheatmap",
  "RColorBrewer",
  "gprofiler2"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
}


library(ggplot2)
library(dplyr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(gprofiler2)


############################################################
# 2. HELPER FUNCTION
#
# g:Profiler can return list-columns.
#
# write.csv() cannot directly write list-columns.
#
# This converts list-columns to semicolon-separated text.
############################################################

flatten_for_csv <- function(df) {
  
  df <- as.data.frame(df)
  
  for (i in seq_along(df)) {
    
    if (is.list(df[[i]])) {
      
      df[[i]] <- vapply(
        
        df[[i]],
        
        function(x) {
          
          if (
            is.null(x) ||
            length(x) == 0
          ) {
            
            return("")
            
          }
          
          paste(
            as.character(
              unlist(
                x,
                use.names = FALSE
              )
            ),
            collapse = ";"
          )
          
        },
        
        character(1)
        
      )
    }
  }
  
  return(df)
}


############################################################
# 3. INPUT / OUTPUT FILES
############################################################

input_file <- "/home/andrew/DESeq2_rhythm/DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv"

background_file <- "/home/andrew/DESeq2_rhythm/DESeq2_24h_rhythm_all_genes.csv"

vst_file <- "/home/andrew/DESeq2_rhythm/DESeq2_VST_expression.csv"

output_dir <- "/home/andrew/DESeq2_rhythm/biological_interpretation"


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 4. READ RHYTHMIC RESULTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "Reading rhythmic gene results\n"
)

cat(
  "====================================================\n"
)


if (!file.exists(input_file)) {
  
  stop(
    "\nInput file does not exist:\n",
    input_file
  )
}


rhythmic <- read_csv(
  input_file,
  show_col_types = FALSE
)


############################################################
# IMPORTANT:
#
# First column contains gene IDs.
#
# Explicitly rename column 1 to "gene".
#
# No automatic detection.
############################################################

colnames(rhythmic)[1] <- "gene"


rhythmic$gene <- as.character(
  rhythmic$gene
)


cat(
  "\nGene IDs are being read from column 1.\n"
)


cat(
  "Number of input rows:",
  nrow(rhythmic),
  "\n"
)


cat(
  "\nInput columns:\n"
)

print(
  colnames(rhythmic)
)


cat(
  "\nFirst five gene IDs:\n"
)

print(
  head(
    rhythmic$gene,
    5
  )
)


############################################################
# 5. CHECK REQUIRED COLUMNS
############################################################

required_columns <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "beta_cos",
  "beta_sin",
  "amplitude",
  "peak_ZT"
)


missing_columns <- setdiff(
  required_columns,
  colnames(rhythmic)
)


if (
  length(missing_columns) > 0
) {
  
  stop(
    
    "\nThe following required columns are missing:\n",
    
    paste(
      missing_columns,
      collapse = ", "
    ),
    
    "\n\nColumns found:\n",
    
    paste(
      colnames(rhythmic),
      collapse = ", "
    )
    
  )
}


############################################################
# 6. REMOVE BLANK GENE IDS
############################################################

rhythmic <- rhythmic %>%
  
  filter(
    !is.na(gene),
    gene != ""
  )


############################################################
# 7. FILTER RHYTHMIC GENES
#
# FDR < 0.05
############################################################

rhythmic <- rhythmic %>%
  
  filter(
    
    !is.na(padj),
    
    padj < 0.05,
    
    !is.na(peak_ZT),
    
    !is.na(amplitude)
    
  )


cat(
  "\nNumber of rhythmic genes (FDR < 0.05):",
  nrow(rhythmic),
  "\n"
)


if (
  nrow(rhythmic) == 0
) {
  
  stop(
    "\nNo genes remain after filtering FDR < 0.05.\n"
  )
}


############################################################
# 8. BASIC SUMMARY
############################################################

summary_table <- data.frame(
  
  Metric = c(
    "Number of rhythmic genes",
    "Minimum FDR",
    "Median FDR",
    "Maximum amplitude",
    "Median amplitude",
    "Minimum peak ZT",
    "Maximum peak ZT"
  ),
  
  Value = c(
    
    nrow(rhythmic),
    
    min(
      rhythmic$padj,
      na.rm = TRUE
    ),
    
    median(
      rhythmic$padj,
      na.rm = TRUE
    ),
    
    max(
      rhythmic$amplitude,
      na.rm = TRUE
    ),
    
    median(
      rhythmic$amplitude,
      na.rm = TRUE
    ),
    
    min(
      rhythmic$peak_ZT,
      na.rm = TRUE
    ),
    
    max(
      rhythmic$peak_ZT,
      na.rm = TRUE
    )
    
  )
)


write.csv(
  
  summary_table,
  
  file.path(
    output_dir,
    "rhythmic_gene_summary.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# 9. PHASE BINS
############################################################

rhythmic <- rhythmic %>%
  
  mutate(
    
    peak_ZT_wrapped = peak_ZT %% 24,
    
    phase_bin = cut(
      
      peak_ZT_wrapped,
      
      breaks = c(
        0,
        4,
        8,
        12,
        16,
        20,
        24
      ),
      
      labels = c(
        "ZT0-4",
        "ZT4-8",
        "ZT8-12",
        "ZT12-16",
        "ZT16-20",
        "ZT20-24"
      ),
      
      include.lowest = TRUE,
      
      right = FALSE
      
    )
    
  )


############################################################
# 10. PHASE COUNTS
############################################################

phase_counts <- rhythmic %>%
  
  count(
    phase_bin,
    .drop = FALSE
  ) %>%
  
  mutate(
    
    percentage = ifelse(
      sum(n) > 0,
      100 * n / sum(n),
      0
    )
    
  )


write.csv(
  
  phase_counts,
  
  file.path(
    output_dir,
    "rhythmic_genes_by_phase.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# 11. PHASE DISTRIBUTION
############################################################

p_phase <- ggplot(
  
  rhythmic,
  
  aes(
    x = peak_ZT_wrapped
  )
  
) +
  
  geom_histogram(
    
    binwidth = 2,
    
    boundary = 0,
    
    closed = "left"
    
  ) +
  
  scale_x_continuous(
    
    breaks = seq(
      0,
      24,
      4
    ),
    
    limits = c(
      0,
      24
    )
    
  ) +
  
  labs(
    
    title = "Distribution of rhythmic gene phases",
    
    x = "Peak time (ZT)",
    
    y = "Number of rhythmic genes"
    
  ) +
  
  theme_classic(
    base_size = 14
  )


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_gene_phase_distribution.pdf"
  ),
  
  p_phase,
  
  width = 8,
  height = 5
  
)


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_gene_phase_distribution.png"
  ),
  
  p_phase,
  
  width = 8,
  height = 5,
  
  dpi = 300
  
)


############################################################
# 12. PHASE BIN BAR PLOT
############################################################

p_phase_bins <- ggplot(
  
  phase_counts,
  
  aes(
    x = phase_bin,
    y = n
  )
  
) +
  
  geom_col() +
  
  labs(
    
    title = "Rhythmic genes by peak phase",
    
    x = "Peak phase",
    
    y = "Number of rhythmic genes"
    
  ) +
  
  theme_classic(
    base_size = 14
  )


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_genes_by_phase_bin.pdf"
  ),
  
  p_phase_bins,
  
  width = 8,
  height = 5
  
)


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_genes_by_phase_bin.png"
  ),
  
  p_phase_bins,
  
  width = 8,
  height = 5,
  
  dpi = 300
  
)


############################################################
# 13. AMPLITUDE DISTRIBUTION
############################################################

p_amplitude <- ggplot(
  
  rhythmic,
  
  aes(
    x = amplitude
  )
  
) +
  
  geom_histogram(
    bins = 30
  ) +
  
  labs(
    
    title = "Amplitude of rhythmic genes",
    
    x = "Cosinor amplitude",
    
    y = "Number of rhythmic genes"
    
  ) +
  
  theme_classic(
    base_size = 14
  )


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_gene_amplitude_distribution.pdf"
  ),
  
  p_amplitude,
  
  width = 8,
  height = 5
  
)


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_gene_amplitude_distribution.png"
  ),
  
  p_amplitude,
  
  width = 8,
  height = 5,
  
  dpi = 300
  
)


############################################################
# 14. PHASE VS AMPLITUDE
############################################################

p_phase_amplitude <- ggplot(
  
  rhythmic,
  
  aes(
    x = peak_ZT_wrapped,
    y = amplitude
  )
  
) +
  
  geom_point(
    alpha = 0.7
  ) +
  
  scale_x_continuous(
    
    breaks = seq(
      0,
      24,
      4
    ),
    
    limits = c(
      0,
      24
    )
    
  ) +
  
  labs(
    
    title = "Rhythmic gene phase and amplitude",
    
    x = "Peak time (ZT)",
    
    y = "Amplitude"
    
  ) +
  
  theme_classic(
    base_size = 14
  )


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_phase_vs_amplitude.pdf"
  ),
  
  p_phase_amplitude,
  
  width = 8,
  height = 5
  
)


ggsave(
  
  file.path(
    output_dir,
    "rhythmic_phase_vs_amplitude.png"
  ),
  
  p_phase_amplitude,
  
  width = 8,
  height = 5,
  
  dpi = 300
  
)


############################################################
# 15. VST PHASE-ORDERED HEATMAP
############################################################

if (
  file.exists(vst_file)
) {
  
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "Creating phase-ordered VST heatmap\n"
  )
  
  cat(
    "====================================================\n"
  )
  
  
  vst <- read.csv(
    
    vst_file,
    
    row.names = 1,
    
    check.names = FALSE
    
  )
  
  
  rownames(vst) <- as.character(
    rownames(vst)
  )
  
  
  common_genes <- intersect(
    
    rhythmic$gene,
    
    rownames(vst)
    
  )
  
  
  cat(
    
    "\nRhythmic genes present in VST matrix:",
    
    length(common_genes),
    
    "\n"
    
  )
  
  
  if (
    length(common_genes) > 1
  ) {
    
    
    rhythmic_ordered <- rhythmic %>%
      
      filter(
        gene %in% common_genes
      ) %>%
      
      arrange(
        peak_ZT_wrapped
      )
    
    
    gene_order <- rhythmic_ordered$gene
    
    
    heatmap_matrix <- as.matrix(
      
      vst[
        gene_order,
        ,
        drop = FALSE
      ]
      
    )
    
    
    ########################################################
    # Row-wise Z-score
    ########################################################
    
    heatmap_matrix <- t(
      
      scale(
        t(
          heatmap_matrix
        )
      )
      
    )
    
    
    ########################################################
    # Sample annotation
    ########################################################
    
    sample_names <- colnames(
      heatmap_matrix
    )
    
    
    expected_ZT <- rep(
      
      c(
        0,
        4,
        8,
        12,
        16,
        20
      ),
      
      each = 4
      
    )
    
    
    if (
      length(expected_ZT) ==
      length(sample_names)
    ) {
      
      
      annotation_col <- data.frame(
        
        ZT = factor(
          
          expected_ZT,
          
          levels = c(
            0,
            4,
            8,
            12,
            16,
            20
          )
          
        )
        
      )
      
      
      rownames(annotation_col) <-
        sample_names
      
      
    } else {
      
      
      warning(
        
        "\nVST matrix contains ",
        
        length(sample_names),
        
        " samples.\n",
        
        "Expected 24 samples.\n",
        
        "Generating heatmap without ZT annotation.\n"
        
      )
      
      
      annotation_col <- NULL
      
    }
    
    
    ########################################################
    # Heatmap PDF
    ########################################################
    
    pdf(
      
      file.path(
        output_dir,
        "all_rhythmic_genes_phase_ordered_heatmap.pdf"
      ),
      
      width = 10,
      height = 16
      
    )
    
    
    pheatmap(
      
      heatmap_matrix,
      
      cluster_rows = FALSE,
      
      cluster_cols = FALSE,
      
      show_rownames = FALSE,
      
      annotation_col = annotation_col,
      
      scale = "none",
      
      main = "FDR < 0.05 rhythmic genes ordered by peak phase"
      
    )
    
    
    dev.off()
    
    
    ########################################################
    # Save matrix
    ########################################################
    
    write.csv(
      
      heatmap_matrix,
      
      file.path(
        output_dir,
        "phase_ordered_rhythmic_gene_VST_matrix.csv"
      )
      
    )
    
    
  } else {
    
    
    message(
      
      "\nNot enough rhythmic genes found in VST matrix.\n",
      
      "Skipping heatmap.\n"
      
    )
    
  }
  
  
} else {
  
  
  message(
    
    "\nVST file not found:\n",
    
    vst_file,
    
    "\nSkipping heatmap.\n"
    
  )
  
}


############################################################
# 16. PREPARE GENE LIST FOR g:PROFILER
############################################################

rhythmic_genes <- unique(
  rhythmic$gene
)


rhythmic_genes <- rhythmic_genes[
  
  !is.na(rhythmic_genes) &
    rhythmic_genes != ""
  
]


cat(
  
  "\nGenes submitted to g:Profiler:",
  
  length(rhythmic_genes),
  
  "\n"
  
)


############################################################
# 17. READ BACKGROUND GENES
#
# IMPORTANT:
#   First column contains gene IDs.
############################################################

background_genes <- NULL


if (
  file.exists(background_file)
) {
  
  
  cat(
    "\nReading background gene list...\n"
  )
  
  
  all_results <- read.csv(
    
    background_file,
    
    check.names = FALSE
    
  )
  
  
  ##########################################################
  # Explicitly use FIRST COLUMN
  ##########################################################
  
  background_genes <- as.character(
    all_results[[1]]
  )
  
  
  background_genes <- unique(
    
    background_genes[
      
      !is.na(background_genes) &
        background_genes != ""
      
    ]
    
  )
  
  
  cat(
    
    "Background genes:",
    
    length(background_genes),
    
    "\n"
    
  )
  
  
} else {
  
  
  message(
    
    "\nBackground file not found:\n",
    
    background_file,
    
    "\n\nRunning g:Profiler without custom background.\n"
    
  )
  
}


############################################################
# 18. RUN OVERALL g:PROFILER
#
# IMPORTANT:
#
# significant = FALSE
#
# This returns all tested terms instead of only terms
# passing FDR < 0.05.
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "Running overall g:Profiler enrichment\n"
)

cat(
  "====================================================\n"
)


if (
  is.null(background_genes)
) {
  
  
  gp <- tryCatch(
    
    
    gost(
      
      query = rhythmic_genes,
      
      organism = "mmusculus",
      
      sources = c(
        "GO:BP",
        "KEGG",
        "REAC"
      ),
      
      significant = FALSE,
      
      correction_method = "fdr",
      
      user_threshold = 0.05,
      
      evcodes = TRUE
      
    ),
    
    
    error = function(e) {
      
      message(
        
        "\ng:Profiler failed:\n",
        
        e$message
        
      )
      
      return(NULL)
      
    }
    
  )
  
  
} else {
  
  
  gp <- tryCatch(
    
    
    gost(
      
      query = rhythmic_genes,
      
      organism = "mmusculus",
      
      sources = c(
        "GO:BP",
        "KEGG",
        "REAC"
      ),
      
      significant = FALSE,
      
      correction_method = "fdr",
      
      user_threshold = 0.05,
      
      custom_bg = background_genes,
      
      evcodes = TRUE
      
    ),
    
    
    error = function(e) {
      
      message(
        
        "\ng:Profiler failed:\n",
        
        e$message
        
      )
      
      return(NULL)
      
    }
    
  )
  
}


############################################################
# 19. PROCESS OVERALL ENRICHMENT
############################################################

if (
  
  !is.null(gp) &&
  
  !is.null(gp$result) &&
  
  nrow(gp$result) > 0
  
) {
  
  
  enrichment_table <- gp$result
  
  
  ##########################################################
  # DIAGNOSTIC: DATABASES RETURNED
  ##########################################################
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "g:Profiler database summary\n"
  )
  
  cat(
    "====================================================\n"
  )
  
  
  print(
    
    table(
      
      enrichment_table$source,
      
      useNA = "ifany"
      
    )
    
  )
  
  
  cat(
    "\nUnique source names returned:\n"
  )
  
  print(
    unique(
      enrichment_table$source
    )
  )
  
  
  ##########################################################
  # Convert list columns
  ##########################################################
  
  enrichment_table_csv <- flatten_for_csv(
    enrichment_table
  )
  
  
  ##########################################################
  # Save ALL enrichment results
  ##########################################################
  
  write.csv(
    
    enrichment_table_csv,
    
    file.path(
      output_dir,
      "gProfiler_overall_enrichment_ALL.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ##########################################################
  # GO
  ##########################################################
  
  go_table <- enrichment_table_csv %>%
    
    filter(
      source == "GO:BP"
    )
  
  
  write.csv(
    
    go_table,
    
    file.path(
      output_dir,
      "GO_Biological_Process_enrichment_all.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  go_significant <- go_table %>%
    
    filter(
      p_value < 0.05
    )
  
  
  write.csv(
    
    go_significant,
    
    file.path(
      output_dir,
      "GO_Biological_Process_enrichment_FDR05.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ##########################################################
  # KEGG
  ##########################################################
  
  kegg_table <- enrichment_table_csv %>%
    
    filter(
      source == "KEGG"
    )
  
  
  cat(
    
    "\nNumber of KEGG terms returned:",
    
    nrow(kegg_table),
    
    "\n"
    
  )
  
  
  if (
    nrow(kegg_table) > 0
  ) {
    
    cat(
      "\nTop KEGG terms before FDR filtering:\n"
    )
    
    
    print(
      
      kegg_table %>%
        
        select(
          term_id,
          term_name,
          p_value,
          intersection_size
        ) %>%
        
        arrange(
          p_value
        ) %>%
        
        head(20)
      
    )
    
  } else {
    
    cat(
      "\nNo KEGG terms were returned by g:Profiler.\n"
    )
    
  }
  
  
  ##########################################################
  # Save ALL KEGG
  ##########################################################
  
  write.csv(
    
    kegg_table,
    
    file.path(
      output_dir,
      "KEGG_enrichment_all.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ##########################################################
  # Significant KEGG
  ##########################################################
  
  kegg_significant <- kegg_table %>%
    
    filter(
      p_value < 0.05
    )
  
  
  write.csv(
    
    kegg_significant,
    
    file.path(
      output_dir,
      "KEGG_enrichment_FDR05.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  cat(
    
    "Number of KEGG terms with FDR < 0.05:",
    
    nrow(kegg_significant),
    
    "\n"
    
  )
  
  
  ##########################################################
  # Reactome
  ##########################################################
  
  reactome_table <- enrichment_table_csv %>%
    
    filter(
      source == "REAC"
    )
  
  
  write.csv(
    
    reactome_table,
    
    file.path(
      output_dir,
      "Reactome_enrichment_all.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  reactome_significant <- reactome_table %>%
    
    filter(
      p_value < 0.05
    )
  
  
  write.csv(
    
    reactome_significant,
    
    file.path(
      output_dir,
      "Reactome_enrichment_FDR05.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ##########################################################
  # Backward-compatible overall file
  #
  # This contains all terms.
  ##########################################################
  
  write.csv(
    
    enrichment_table_csv,
    
    file.path(
      output_dir,
      "gProfiler_overall_enrichment.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ##########################################################
  # OVERALL ENRICHMENT PLOT
  #
  # Plot top 20 terms by adjusted P value.
  ##########################################################
  
  plot_data <- enrichment_table %>%
    
    arrange(
      p_value
    ) %>%
    
    slice_head(
      n = 20
    )
  
  
  if (
    nrow(plot_data) > 0
  ) {
    
    
    plot_data$term_name <- factor(
      
      plot_data$term_name,
      
      levels = rev(
        plot_data$term_name
      )
      
    )
    
    
    p_enrichment <- ggplot(
      
      plot_data,
      
      aes(
        x = -log10(p_value),
        y = term_name
      )
      
    ) +
      
      geom_point(
        
        aes(
          size = intersection_size
        )
        
      ) +
      
      labs(
        
        title = "Top enriched pathways/processes",
        
        x = "-log10(FDR-adjusted P value)",
        
        y = "Pathway / biological process",
        
        size = "Genes"
        
      ) +
      
      theme_classic(
        base_size = 12
      )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "gProfiler_overall_enrichment.pdf"
      ),
      
      p_enrichment,
      
      width = 10,
      height = 8
      
    )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "gProfiler_overall_enrichment.png"
      ),
      
      p_enrichment,
      
      width = 10,
      height = 8,
      
      dpi = 300
      
    )
    
  }
  
  
} else {
  
  
  message(
    
    "\nNo enrichment results were returned by g:Profiler.\n"
    
  )
  
}


############################################################
# 20. PHASE-SPECIFIC ENRICHMENT
#
# Again use significant = FALSE so that ALL terms are
# available for inspection.
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "Running phase-specific enrichment\n"
)

cat(
  "====================================================\n"
)


phase_levels <- levels(
  rhythmic$phase_bin
)


phase_enrichment_summary <- list()


for (
  phase in phase_levels
) {
  
  
  cat(
    "\n----------------------------------------------------\n"
  )
  
  cat(
    "Phase:",
    phase,
    "\n"
  )
  
  
  ##########################################################
  # Genes in phase
  ##########################################################
  
  genes_phase <- rhythmic %>%
    
    filter(
      phase_bin == phase
    ) %>%
    
    pull(
      gene
    )
  
  
  genes_phase <- unique(
    
    genes_phase[
      
      !is.na(genes_phase) &
        genes_phase != ""
      
    ]
    
  )
  
  
  cat(
    
    "Number of genes:",
    
    length(genes_phase),
    
    "\n"
    
  )
  
  
  ##########################################################
  # Minimum number of genes
  ##########################################################
  
  if (
    length(genes_phase) < 5
  ) {
    
    
    message(
      
      "Skipping ",
      
      phase,
      
      ": fewer than 5 genes."
      
    )
    
    
    next
    
  }
  
  
  ##########################################################
  # Run g:Profiler
  ##########################################################
  
  if (
    is.null(background_genes)
  ) {
    
    
    gp_phase <- tryCatch(
      
      
      gost(
        
        query = genes_phase,
        
        organism = "mmusculus",
        
        sources = c(
          "GO:BP",
          "KEGG",
          "REAC"
        ),
        
        significant = FALSE,
        
        correction_method = "fdr",
        
        user_threshold = 0.05,
        
        evcodes = TRUE
        
      ),
      
      
      error = function(e) {
        
        message(
          
          "g:Profiler failed for ",
          
          phase,
          
          ": ",
          
          e$message
          
        )
        
        return(NULL)
        
      }
      
    )
    
    
  } else {
    
    
    gp_phase <- tryCatch(
      
      
      gost(
        
        query = genes_phase,
        
        organism = "mmusculus",
        
        sources = c(
          "GO:BP",
          "KEGG",
          "REAC"
        ),
        
        significant = FALSE,
        
        correction_method = "fdr",
        
        user_threshold = 0.05,
        
        custom_bg = background_genes,
        
        evcodes = TRUE
        
      ),
      
      
      error = function(e) {
        
        message(
          
          "g:Profiler failed for ",
          
          phase,
          
          ": ",
          
          e$message
          
        )
        
        return(NULL)
        
      }
      
    )
    
  }
  
  
  ##########################################################
  # Save results
  ##########################################################
  
  if (
    
    !is.null(gp_phase) &&
    
    !is.null(gp_phase$result) &&
    
    nrow(gp_phase$result) > 0
    
  ) {
    
    
    phase_table <- gp_phase$result
    
    
    phase_table$Phase <- phase
    
    
    ########################################################
    # Diagnostic
    ########################################################
    
    cat(
      
      "Databases returned:\n"
      
    )
    
    print(
      
      table(
        
        phase_table$source,
        
        useNA = "ifany"
        
      )
      
    )
    
    
    ########################################################
    # Flatten list columns
    ########################################################
    
    phase_table_csv <- flatten_for_csv(
      phase_table
    )
    
    
    ########################################################
    # Save all phase enrichment
    ########################################################
    
    write.csv(
      
      phase_table_csv,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_enrichment_all.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    ########################################################
    # GO
    ########################################################
    
    phase_go <- phase_table_csv %>%
      
      filter(
        source == "GO:BP"
      )
    
    
    phase_go_fdr <- phase_go %>%
      
      filter(
        p_value < 0.05
      )
    
    
    write.csv(
      
      phase_go,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_GO_BP_all.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    write.csv(
      
      phase_go_fdr,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_GO_BP_FDR05.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    ########################################################
    # KEGG
    ########################################################
    
    phase_kegg <- phase_table_csv %>%
      
      filter(
        source == "KEGG"
      )
    
    
    phase_kegg_fdr <- phase_kegg %>%
      
      filter(
        p_value < 0.05
      )
    
    
    write.csv(
      
      phase_kegg,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_KEGG_all.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    write.csv(
      
      phase_kegg_fdr,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_KEGG_FDR05.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    cat(
      
      "KEGG terms:",
      
      nrow(phase_kegg),
      
      "\n"
      
    )
    
    
    cat(
      
      "KEGG terms FDR < 0.05:",
      
      nrow(phase_kegg_fdr),
      
      "\n"
      
    )
    
    
    ########################################################
    # Reactome
    ########################################################
    
    phase_reactome <- phase_table_csv %>%
      
      filter(
        source == "REAC"
      )
    
    
    phase_reactome_fdr <- phase_reactome %>%
      
      filter(
        p_value < 0.05
      )
    
    
    write.csv(
      
      phase_reactome,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_Reactome_all.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    write.csv(
      
      phase_reactome_fdr,
      
      file.path(
        output_dir,
        paste0(
          "gProfiler_",
          phase,
          "_Reactome_FDR05.csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
    
    ########################################################
    # Store combined phase table
    ########################################################
    
    phase_enrichment_summary[[phase]] <-
      phase_table_csv
    
    
    ########################################################
    # Plot top 15 terms
    ########################################################
    
    plot_phase <- phase_table %>%
      
      arrange(
        p_value
      ) %>%
      
      slice_head(
        n = 15
      )
    
    
    if (
      nrow(plot_phase) > 0
    ) {
      
      
      plot_phase$term_name <- factor(
        
        plot_phase$term_name,
        
        levels = rev(
          plot_phase$term_name
        )
        
      )
      
      
      p_phase_enrichment <- ggplot(
        
        plot_phase,
        
        aes(
          x = -log10(p_value),
          y = term_name
        )
        
      ) +
        
        geom_point(
          
          aes(
            size = intersection_size
          )
          
        ) +
        
        labs(
          
          title = paste(
            "Top enrichment:",
            phase
          ),
          
          x = "-log10(FDR-adjusted P value)",
          
          y = "Pathway / biological process",
          
          size = "Genes"
          
        ) +
        
        theme_classic(
          base_size = 12
        )
      
      
      ggsave(
        
        file.path(
          output_dir,
          paste0(
            "gProfiler_",
            phase,
            "_enrichment.pdf"
          )
        ),
        
        p_phase_enrichment,
        
        width = 10,
        height = 7
        
      )
      
      
      ggsave(
        
        file.path(
          output_dir,
          paste0(
            "gProfiler_",
            phase,
            "_enrichment.png"
          )
        ),
        
        p_phase_enrichment,
        
        width = 10,
        height = 7,
        
        dpi = 300
        
      )
      
    }
    
    
  } else {
    
    
    message(
      
      "No enrichment results returned for ",
      
      phase
      
    )
    
  }
  
}


############################################################
# 21. COMBINED PHASE ENRICHMENT
############################################################

if (
  length(phase_enrichment_summary) > 0
) {
  
  
  combined_phase_enrichment <- bind_rows(
    phase_enrichment_summary
  )
  
  
  write.csv(
    
    combined_phase_enrichment,
    
    file.path(
      output_dir,
      "gProfiler_phase_enrichment_combined.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
} else {
  
  
  message(
    
    "\nNo phase-specific enrichment results were generated.\n"
    
  )
  
}


############################################################
# 22. FINAL RHYTHMIC GENE TABLE
############################################################

rhythmic_final <- rhythmic %>%
  
  arrange(
    peak_ZT_wrapped,
    padj
  )


write.csv(
  
  rhythmic_final,
  
  file.path(
    output_dir,
    "rhythmic_genes_final_phase_ordered.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# 23. TOP 10 GENES PER PHASE
############################################################

top_by_phase <- rhythmic %>%
  
  group_by(
    phase_bin
  ) %>%
  
  arrange(
    padj
  ) %>%
  
  slice_head(
    n = 10
  ) %>%
  
  ungroup()


write.csv(
  
  top_by_phase,
  
  file.path(
    output_dir,
    "top_10_rhythmic_genes_per_phase.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# 24. FINAL OUTPUT SUMMARY
############################################################

cat(
  "\n\n"
)

cat(
  "====================================================\n"
)

cat(
  "ANALYSIS COMPLETE\n"
)

cat(
  "====================================================\n"
)


cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)


cat(
  "\nNumber of rhythmic genes:",
  nrow(rhythmic),
  "\n"
)


cat(
  "\nMain enrichment outputs:\n"
)

cat(
  "  gProfiler_overall_enrichment_ALL.csv\n"
)

cat(
  "  gProfiler_overall_enrichment.csv\n"
)

cat(
  "  GO_Biological_Process_enrichment_all.csv\n"
)

cat(
  "  GO_Biological_Process_enrichment_FDR05.csv\n"
)

cat(
  "  KEGG_enrichment_all.csv\n"
)

cat(
  "  KEGG_enrichment_FDR05.csv\n"
)

cat(
  "  Reactome_enrichment_all.csv\n"
)

cat(
  "  Reactome_enrichment_FDR05.csv\n"
)


cat(
  "\nPhase-specific enrichment files include:\n"
)

cat(
  "  *_GO_BP_all.csv\n"
)

cat(
  "  *_GO_BP_FDR05.csv\n"
)

cat(
  "  *_KEGG_all.csv\n"
)

cat(
  "  *_KEGG_FDR05.csv\n"
)

cat(
  "  *_Reactome_all.csv\n"
)

cat(
  "  *_Reactome_FDR05.csv\n"
)


cat(
  "\n====================================================\n"
)

cat(
  "Finished.\n"
)

cat(
  "====================================================\n"
)
