############################################################
# DESeq2 24-HOUR RHYTHMIC GENE EXPRESSION ANALYSIS
#
# Experiment:
#   6 timepoints:
#       ZT0, ZT4, ZT8, ZT12, ZT16, ZT20
#
#   4 biological replicates per timepoint
#
# Model:
#   Full    = ~ cos24 + sin24
#   Reduced = ~ 1
#
# Statistical test:
#   DESeq2 likelihood-ratio test (LRT)
#
# Cosinor plots:
#   ONLY the TOP_N genes with FDR < 0.05
#   ranked by lowest adjusted p-value
#
############################################################


############################
# 1. Load packages
############################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})


############################
# 2. Input/output settings
############################

count_file <- "/home/andrew/NMR_featureCounts_clean.txt"

outdir <- "/home/andrew/DESeq2_rhythm"

# Maximum number of significant genes to plot
TOP_N <- 50


dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################
# 3. Initialize objects
############################

# These are initialized so that the final summary
# cannot fail if no genes pass FDR < 0.05.

rhythmic_results <- data.frame()

top_genes <- character(0)


############################
# 4. Read featureCounts data
############################

cat("\nReading featureCounts data...\n")

counts_raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat(
  "Count table dimensions: ",
  nrow(counts_raw),
  " genes x ",
  ncol(counts_raw),
  " columns\n",
  sep = ""
)


############################
# 5. Check featureCounts columns
############################

if (ncol(counts_raw) < 7) {
  
  stop(
    "The featureCounts file has fewer than 7 columns. ",
    "Expected Geneid + 5 annotation columns + sample counts."
  )
  
}


if (!"Geneid" %in% colnames(counts_raw)) {
  
  stop(
    "Could not find a 'Geneid' column in the featureCounts file."
  )
  
}


############################
# 6. Extract count matrix
############################

# featureCounts output:
#
# Column 1 = Geneid
# Columns 2-6 = annotation information
# Columns 7 onwards = sample counts

count_matrix <- counts_raw[
  ,
  7:ncol(counts_raw),
  drop = FALSE
]


rownames(count_matrix) <- counts_raw$Geneid


# Check for duplicated gene IDs

if (anyDuplicated(rownames(count_matrix)) > 0) {
  
  duplicated_genes <- unique(
    rownames(count_matrix)[
      duplicated(rownames(count_matrix))
    ]
  )
  
  stop(
    "Duplicated Geneid values detected. ",
    "Please resolve duplicated gene IDs before running DESeq2.\n",
    "Examples: ",
    paste(
      head(duplicated_genes, 10),
      collapse = ", "
    )
  )
  
}


count_matrix <- as.matrix(
  count_matrix
)


# Convert safely to numeric

mode(count_matrix) <- "numeric"


# Check for NA values

if (anyNA(count_matrix)) {
  
  stop(
    "NA values detected in the count matrix."
  )
  
}


# Counts must be non-negative

if (any(count_matrix < 0)) {
  
  stop(
    "Negative count values detected."
  )
  
}


# Round if necessary

if (any(count_matrix != round(count_matrix))) {
  
  warning(
    "Non-integer counts detected. ",
    "Counts will be rounded to integers."
  )
  
}


storage.mode(count_matrix) <- "integer"


############################
# 7. Sample metadata
############################

sample_info <- data.frame(
  
  sample = c(
    
    # ZT0
    "NMR2-COL1-HYP-13",
    "NMR2-COL1-HYP-14",
    "NMR2-COL1-HYP-15",
    "NMR2-COL1-HYP-16",
    
    # ZT4
    "NMR2-COL1-HYP-01",
    "NMR2-COL1-HYP-02",
    "NMR2-COL1-HYP-03",
    "NMR2-COL1-HYP-04",
    
    # ZT8
    "NMR2-COL1-HYP-05",
    "NMR2-COL1-HYP-06",
    "NMR2-COL1-HYP-07",
    "NMR2-COL1-HYP-08",
    
    # ZT12
    "NMR2-COL1-HYP-09",
    "NMR2-COL1-HYP-10",
    "NMR2-COL1-HYP-11",
    "NMR2-COL1-HYP-12",
    
    # ZT16
    "NMR2-COL1-HYP-17",
    "NMR2-COL1-HYP-18",
    "NMR2-COL1-HYP-19",
    "NMR2-COL1-HYP-20",
    
    # ZT20
    "NMR2-COL1-HYP-21",
    "NMR2-COL1-HYP-22",
    "NMR2-COL1-HYP-23",
    "NMR2-COL1-HYP-24"
  ),
  
  ZT = c(
    0, 0, 0, 0,
    4, 4, 4, 4,
    8, 8, 8, 8,
    12, 12, 12, 12,
    16, 16, 16, 16,
    20, 20, 20, 20
  ),
  
  replicate = c(
    1, 2, 3, 4,
    1, 2, 3, 4,
    1, 2, 3, 4,
    1, 2, 3, 4,
    1, 2, 3, 4,
    1, 2, 3, 4
  ),
  
  stringsAsFactors = FALSE
)


rownames(sample_info) <- sample_info$sample


############################
# 8. Check metadata
############################

if (nrow(sample_info) != 24) {
  
  stop(
    "Expected 24 samples in metadata."
  )
  
}


if (anyDuplicated(sample_info$sample)) {
  
  stop(
    "Duplicated sample names detected in metadata."
  )
  
}


############################
# 9. Check sample names
############################

cat("\nChecking sample names...\n")


missing_from_counts <- setdiff(
  rownames(sample_info),
  colnames(count_matrix)
)


missing_from_metadata <- setdiff(
  colnames(count_matrix),
  rownames(sample_info)
)


if (length(missing_from_counts) > 0) {
  
  cat(
    "\nSamples in metadata but missing from count matrix:\n"
  )
  
  print(missing_from_counts)
  
  stop(
    "Sample names do not match."
  )
  
}


if (length(missing_from_metadata) > 0) {
  
  cat(
    "\nSamples in count matrix but missing from metadata:\n"
  )
  
  print(missing_from_metadata)
  
  stop(
    "Sample names do not match."
  )
  
}


# Reorder count matrix to match metadata exactly

count_matrix <- count_matrix[
  ,
  rownames(sample_info),
  drop = FALSE
]


stopifnot(
  identical(
    colnames(count_matrix),
    rownames(sample_info)
  )
)


cat(
  "Sample names match successfully.\n"
)


############################
# 10. Create harmonic variables
############################

omega <- 2 * pi / 24


sample_info$cos24 <- cos(
  omega * sample_info$ZT
)


sample_info$sin24 <- sin(
  omega * sample_info$ZT
)


############################
# 11. Save metadata
############################

write.table(
  
  sample_info,
  
  file = file.path(
    outdir,
    "sample_metadata.tsv"
  ),
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = TRUE,
  
  col.names = NA
)


############################
# 12. Filter low-count genes
############################

cat(
  "\nFiltering low-count genes...\n"
)


# Keep genes with:
# >=10 reads in at least 4 samples

keep <- rowSums(
  count_matrix >= 10
) >= 4


cat(
  "Genes before filtering: ",
  nrow(count_matrix),
  "\n",
  sep = ""
)


cat(
  "Genes after filtering: ",
  sum(keep),
  "\n",
  sep = ""
)


if (sum(keep) == 0) {
  
  stop(
    "No genes passed the filtering threshold."
  )
  
}


count_matrix_filtered <- count_matrix[
  keep,
  ,
  drop = FALSE
]


############################
# 13. Create DESeqDataSet
############################

dds <- DESeqDataSetFromMatrix(
  
  countData = count_matrix_filtered,
  
  colData = sample_info,
  
  design = ~ cos24 + sin24
)


############################
# 14. Run DESeq2 LRT
############################

cat(
  "\nRunning DESeq2 likelihood-ratio test...\n"
)


dds <- DESeq(
  
  dds,
  
  test = "LRT",
  
  reduced = ~ 1
)


############################
# 15. Display model
############################

cat(
  "\nDESeq2 design:\n"
)

print(
  design(dds)
)


cat(
  "\nDESeq2 coefficient names:\n"
)

print(
  resultsNames(dds)
)


############################
# 16. Identify coefficients
############################

coef_names <- resultsNames(dds)


if (!"cos24" %in% coef_names) {
  
  stop(
    "cos24 coefficient not found in DESeq2 model."
  )
  
}


if (!"sin24" %in% coef_names) {
  
  stop(
    "sin24 coefficient not found in DESeq2 model."
  )
  
}


# Find the intercept automatically

intercept_name <- coef_names[
  grepl(
    "Intercept",
    coef_names,
    ignore.case = TRUE
  )
]


if (length(intercept_name) != 1) {
  
  stop(
    "Could not uniquely identify the DESeq2 intercept coefficient."
  )
  
}


intercept_name <- intercept_name[1]


cat(
  "\nIntercept coefficient: ",
  intercept_name,
  "\n",
  sep = ""
)


############################
# 17. Extract LRT results
############################

res <- results(
  dds
)


# Convert to data frame

res_df <- as.data.frame(
  res
)


# Add gene IDs

res_df$gene_id <- rownames(
  res_df
)


# Rank by adjusted p-value

res_df <- res_df[
  order(
    res_df$padj,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


############################
# 18. Save ALL results
############################

write.csv(
  
  res_df,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_all_genes.csv"
  ),
  
  row.names = TRUE
)


############################
# 19. Select significant genes
############################

rhythmic <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05,
  ,
  drop = FALSE
]


cat(
  "\n============================================\n"
)


cat(
  "FDR < 0.05 rhythmic genes: ",
  nrow(rhythmic),
  "\n",
  sep = ""
)


cat(
  "============================================\n"
)


############################
# 20. Save FDR < 0.05 results
############################

write.csv(
  
  rhythmic,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythmic_genes_FDR05.csv"
  ),
  
  row.names = TRUE
)


############################
# 21. Extract DESeq2 coefficients
############################

coef_matrix <- coef(
  dds
)


cat(
  "\nCoefficient matrix columns:\n"
)

print(
  colnames(coef_matrix)
)


beta_cos <- coef_matrix[
  ,
  "cos24"
]


beta_sin <- coef_matrix[
  ,
  "sin24"
]


gene_intercept_all <- coef_matrix[
  ,
  intercept_name
]


############################
# 22. Calculate amplitude
############################

# For:
#
# y = intercept +
#     beta_cos*cos(wt) +
#     beta_sin*sin(wt)
#
# amplitude = sqrt(
#   beta_cos^2 + beta_sin^2
# )

amplitude <- sqrt(
  
  beta_cos^2 +
    
    beta_sin^2
)


############################
# 23. Calculate peak phase
############################

# Maximum occurs when:
#
# phase = atan2(beta_sin, beta_cos)
#
# converted to hours.

phase_radians <- atan2(
  
  beta_sin,
  
  beta_cos
)


phase_hours <- (
  
  phase_radians %% (2 * pi)
  
) * 24 / (2 * pi)


############################
# 24. Create complete rhythm table
############################

rhythm_results <- res_df


rhythm_results$beta_cos <- beta_cos[
  match(
    rhythm_results$gene_id,
    names(beta_cos)
  )
]


rhythm_results$beta_sin <- beta_sin[
  match(
    rhythm_results$gene_id,
    names(beta_sin)
  )
]


rhythm_results$amplitude <- amplitude[
  match(
    rhythm_results$gene_id,
    names(amplitude)
  )
]


rhythm_results$peak_ZT <- phase_hours[
  match(
    rhythm_results$gene_id,
    names(phase_hours)
  )
]


rhythm_results$intercept <- gene_intercept_all[
  match(
    rhythm_results$gene_id,
    names(gene_intercept_all)
  )
]


# Keep adjusted-p-value ranking

rhythm_results <- rhythm_results[
  order(
    rhythm_results$padj,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


############################
# 25. Save complete rhythm table
############################

write.csv(
  
  rhythm_results,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_with_phase_amplitude.csv"
  ),
  
  row.names = TRUE
)


############################
# 26. Create significant rhythm table
############################

rhythmic_results <- rhythm_results[
  !is.na(rhythm_results$padj) &
    rhythmic_results$padj < 0.05,
  ,
  drop = FALSE
]


rhythmic_results <- rhythmic_results[
  order(
    rhythmic_results$padj,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


############################
# 27. Save significant rhythm table
############################

write.csv(
  
  rhythmic_results,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv"
  ),
  
  row.names = TRUE
)


############################
# 28. Normalized counts
############################

norm_counts <- counts(
  
  dds,
  
  normalized = TRUE
)


write.csv(
  
  as.data.frame(norm_counts),
  
  file = file.path(
    outdir,
    "DESeq2_normalized_counts.csv"
  ),
  
  row.names = TRUE
)


############################
# 29. VST
############################

cat(
  "\nRunning variance stabilizing transformation...\n"
)


vsd <- vst(
  
  dds,
  
  blind = FALSE
)


vsd_matrix <- assay(
  vsd
)


write.csv(
  
  as.data.frame(vsd_matrix),
  
  file = file.path(
    outdir,
    "DESeq2_VST_expression.csv"
  ),
  
  row.names = TRUE
)


############################
# 30. PCA
############################

pca_data <- plotPCA(
  
  vsd,
  
  intgroup = "ZT",
  
  returnData = TRUE
)


percentVar <- attr(
  
  pca_data,
  
  "percentVar"
)


pca_plot <- ggplot(
  
  pca_data,
  
  aes(
    x = PC1,
    y = PC2,
    label = name,
    color = factor(ZT)
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    vjust = -1,
    size = 3
  ) +
  
  xlab(
    paste0(
      "PC1: ",
      round(
        100 * percentVar[1]
      ),
      "%"
    )
  ) +
  
  ylab(
    paste0(
      "PC2: ",
      round(
        100 * percentVar[2]
      ),
      "%"
    )
  ) +
  
  theme_bw() +
  
  labs(
    color = "ZT"
  )


ggsave(
  
  file.path(
    outdir,
    "PCA_ZT.pdf"
  ),
  
  pca_plot,
  
  width = 8,
  
  height = 6
)


############################
# 31. Sample correlation heatmap
############################

sample_cor <- cor(
  
  vsd_matrix,
  
  method = "pearson"
)


pdf(
  
  file.path(
    outdir,
    "sample_correlation_heatmap.pdf"
  ),
  
  width = 10,
  
  height = 10
)


pheatmap(
  
  sample_cor,
  
  annotation_col =
    data.frame(
      ZT = factor(
        sample_info$ZT
      ),
      row.names = rownames(
        sample_info
      )
    ),
  
  annotation_row =
    data.frame(
      ZT = factor(
        sample_info$ZT
      ),
      row.names = rownames(
        sample_info
      )
    ),
  
  main = "Sample correlation"
)


dev.off()


############################
# 32. MA plot
############################

pdf(
  
  file.path(
    outdir,
    "DESeq2_rhythm_MAplot.pdf"
  ),
  
  width = 8,
  
  height = 6
)


plotMA(
  
  res,
  
  alpha = 0.05
)


dev.off()


############################
# 33. Top 20 rhythmic gene heatmap
############################

if (nrow(rhythmic_results) > 0) {
  
  top_heatmap_genes <- rhythmic_results$gene_id[
    1:min(
      20,
      nrow(rhythmic_results)
    )
  ]
  
  
  top_expression <- vsd_matrix[
    
    top_heatmap_genes,
    
    ,
    
    drop = FALSE
  ]
  
  
  top_scaled <- t(
    
    scale(
      t(top_expression)
    )
  )
  
  
  annotation_col <- data.frame(
    
    ZT = factor(
      
      sample_info$ZT,
      
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
    rownames(sample_info)
  
  
  pdf(
    
    file.path(
      outdir,
      "top_20_rhythmic_genes_heatmap.pdf"
    ),
    
    width = 10,
    
    height = 10
  )
  
  
  pheatmap(
    
    top_scaled,
    
    annotation_col = annotation_col,
    
    cluster_cols = FALSE,
    
    scale = "none",
    
    main = "Top 20 rhythmic genes"
  )
  
  
  dev.off()
  
} else {
  
  cat(
    "\nNo FDR < 0.05 genes available for heatmap.\n"
  )
  
}


############################################################
# 34. COSINOR PLOTS
#
# Only the TOP_N genes from the FDR < 0.05 table
# are plotted.
#
# Genes are ranked by lowest adjusted p-value.
############################################################

cat(
  "\n============================================\n"
)

cat(
  "Preparing top rhythmic genes for cosinor plots\n"
)

cat(
  "============================================\n"
)


if (nrow(rhythmic_results) == 0) {
  
  cat(
    "\nNo genes passed FDR < 0.05.\n"
  )
  
  cat(
    "No cosinor plots will be generated.\n"
  )
  
} else {
  
  
  ##########################################################
  # 34.1 Select top genes
  ##########################################################
  
  top_gene_count <- min(
    
    TOP_N,
    
    nrow(rhythmic_results)
  )
  
  
  top_genes <- rhythmic_results$gene_id[
    1:top_gene_count
  ]
  
  
  top_results <- rhythmic_results[
    
    1:top_gene_count,
    
    ,
    
    drop = FALSE
  ]
  
  
  ##########################################################
  # 34.2 Save top-gene table
  ##########################################################
  
  write.csv(
    
    top_results,
    
    file = file.path(
      
      outdir,
      
      paste0(
        "top_",
        TOP_N,
        "_rhythmic_genes_by_FDR.csv"
      )
    ),
    
    row.names = TRUE
  )
  
  
  ##########################################################
  # 34.3 Cosinor output directory
  ##########################################################
  
  cosinor_dir <- file.path(
    
    outdir,
    
    paste0(
      "top_",
      TOP_N,
      "_cosinor_plots"
    )
  )
  
  
  dir.create(
    
    cosinor_dir,
    
    showWarnings = FALSE,
    
    recursive = TRUE
  )
  
  
  cat(
    "\nGenerating cosinor plots for TOP ",
    length(top_genes),
    " genes.\n",
    sep = ""
  )
  
  
  ##########################################################
  # 34.4 Combined PDF
  ##########################################################
  
  combined_pdf <- file.path(
    
    outdir,
    
    paste0(
      "top_",
      TOP_N,
      "_rhythmic_genes_cosinor_plots.pdf"
    )
  )
  
  
  pdf(
    
    combined_pdf,
    
    width = 8,
    
    height = 6
  )
  
  
  ##########################################################
  # 34.5 Loop through top genes
  ##########################################################
  
  for (gene in top_genes) {
    
    
    cat(
      "Plotting: ",
      gene,
      "\n",
      sep = ""
    )
    
    
    ########################################################
    # Observed normalized counts
    ########################################################
    
    gene_counts <- norm_counts[
      gene,
      ,
      drop = TRUE
    ]
    
    
    plot_data <- data.frame(
      
      sample = rownames(sample_info),
      
      ZT = sample_info$ZT,
      
      replicate = sample_info$replicate,
      
      expression = as.numeric(
        gene_counts
      )
    )
    
    
    ########################################################
    # Mean expression
    ########################################################
    
    mean_data <- aggregate(
      
      expression ~ ZT,
      
      data = plot_data,
      
      FUN = mean
    )
    
    
    ########################################################
    # Standard deviation
    ########################################################
    
    sd_data <- aggregate(
      
      expression ~ ZT,
      
      data = plot_data,
      
      FUN = sd
    )
    
    
    names(sd_data)[2] <- "SD"
    
    
    summary_data <- merge(
      
      mean_data,
      
      sd_data,
      
      by = "ZT"
    )
    
    
    ########################################################
    # Harmonic coefficients
    ########################################################
    
    gene_beta_cos <- beta_cos[
      gene
    ]
    
    
    gene_beta_sin <- beta_sin[
      gene
    ]
    
    
    gene_intercept <- gene_intercept_all[
      gene
    ]
    
    
    ########################################################
    # Check coefficients
    ########################################################
    
    if (
      !is.finite(gene_beta_cos) ||
      !is.finite(gene_beta_sin) ||
      !is.finite(gene_intercept)
    ) {
      
      warning(
        "Non-finite coefficient for ",
        gene,
        ". Skipping."
      )
      
      next
    }
    
    
    ########################################################
    # Smooth 24-hour prediction
    ########################################################
    
    prediction_ZT <- seq(
      
      0,
      
      24,
      
      length.out = 500
    )
    
    
    prediction_cos <- cos(
      
      omega *
        
        prediction_ZT
    )
    
    
    prediction_sin <- sin(
      
      omega *
        
        prediction_ZT
    )
    
    
    ########################################################
    # Linear predictor on log2 scale
    ########################################################
    
    prediction_log2 <- (
      
      gene_intercept +
        
        gene_beta_cos *
        prediction_cos +
        
        gene_beta_sin *
        prediction_sin
    )
    
    
    ########################################################
    # Convert to normalized-expression scale
    ########################################################
    
    prediction_fit <- 2^(
      prediction_log2
    )
    
    
    prediction_data <- data.frame(
      
      ZT = prediction_ZT,
      
      fit = prediction_fit
    )
    
    
    ########################################################
    # Peak phase
    ########################################################
    
    gene_peak <- phase_hours[
      gene
    ]
    
    
    gene_amplitude <- amplitude[
      gene
    ]
    
    
    gene_padj <- rhythmic_results[
      rhythmic_results$gene_id == gene,
      "padj"
    ][1]
    
    
    ########################################################
    # Plot
    ########################################################
    
    p <- ggplot(
      
      plot_data,
      
      aes(
        x = ZT,
        y = expression
      )
      
    ) +
      
      
      ######################################################
    # Individual replicates
    ######################################################
    
    geom_jitter(
      
      width = 0.35,
      
      height = 0,
      
      size = 2.4,
      
      alpha = 0.65
    ) +
      
      
      ######################################################
    # Mean +/- SD
    ######################################################
    
    geom_errorbar(
      
      data = summary_data,
      
      aes(
        
        x = ZT,
        
        ymin = pmax(
          0,
          expression - SD
        ),
        
        ymax =
          expression + SD
      ),
      
      width = 0.4,
      
      linewidth = 0.7
    ) +
      
      
      ######################################################
    # Mean
    ######################################################
    
    geom_point(
      
      data = summary_data,
      
      aes(
        
        x = ZT,
        
        y = expression
      ),
      
      size = 4
    ) +
      
      
      ######################################################
    # Cosinor fitted curve
    ######################################################
    
    geom_line(
      
      data = prediction_data,
      
      aes(
        
        x = ZT,
        
        y = fit
      ),
      
      linewidth = 1.2
    ) +
      
      
      ######################################################
    # Predicted peak
    ######################################################
    
    geom_vline(
      
      xintercept = gene_peak,
      
      linetype = "dashed",
      
      linewidth = 0.7
    ) +
      
      
      ######################################################
    # X axis
    ######################################################
    
    scale_x_continuous(
      
      breaks = c(
        0,
        4,
        8,
        12,
        16,
        20,
        24
      ),
      
      limits = c(
        0,
        24
      )
    ) +
      
      
      ######################################################
    # Y axis
    ######################################################
    
    expand_limits(
      y = 0
    ) +
      
      
      ######################################################
    # Labels
    ######################################################
    
    labs(
      
      title = gene,
      
      subtitle = paste0(
        
        "FDR = ",
        
        format.pval(
          
          gene_padj,
          
          digits = 3,
          
          eps = 1e-10
        ),
        
        "   |   Amplitude = ",
        
        round(
          gene_amplitude,
          3
        ),
        
        " log2 units",
        
        "   |   Peak ZT = ",
        
        round(
          gene_peak,
          2
        )
      ),
      
      x = "Zeitgeber Time (ZT)",
      
      y = "DESeq2-normalized expression"
    ) +
      
      
      ######################################################
    # Theme
    ######################################################
    
    theme_bw(
      
      base_size = 14
    ) +
      
      theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 16
          ),
        
        plot.subtitle =
          element_text(
            size = 11
          ),
        
        axis.title =
          element_text(
            face = "bold"
          )
      )
    
    
    ########################################################
    # Add plot to combined PDF
    ########################################################
    
    print(p)
    
    
    ########################################################
    # Make safe filename
    ########################################################
    
    safe_gene_name <- gsub(
      
      "[^A-Za-z0-9_.-]",
      
      "_",
      
      gene
    )
    
    
    ########################################################
    # Individual PDF
    ########################################################
    
    ggsave(
      
      filename = file.path(
        
        cosinor_dir,
        
        paste0(
          safe_gene_name,
          "_cosinor.pdf"
        )
      ),
      
      plot = p,
      
      width = 8,
      
      height = 6
    )
    
    
    ########################################################
    # Individual PNG
    ########################################################
    
    ggsave(
      
      filename = file.path(
        
        cosinor_dir,
        
        paste0(
          safe_gene_name,
          "_cosinor.png"
        )
      ),
      
      plot = p,
      
      width = 8,
      
      height = 6,
      
      dpi = 300
    )
    
  }
  
  
  ##########################################################
  # Close combined PDF
  ##########################################################
  
  dev.off()
  
  
  ##########################################################
  # Plotting summary
  ##########################################################
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "COSINOR PLOTTING COMPLETE\n"
  )
  
  cat(
    "============================================\n"
  )
  
  cat(
    "FDR < 0.05 genes available: ",
    nrow(rhythmic_results),
    "\n",
    sep = ""
  )
  
  cat(
    "Genes plotted: ",
    length(top_genes),
    "\n",
    sep = ""
  )
  
  cat(
    "\nTop-gene table:\n",
    file.path(
      outdir,
      paste0(
        "top_",
        TOP_N,
        "_rhythmic_genes_by_FDR.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "\nCombined cosinor PDF:\n",
    combined_pdf,
    "\n"
  )
  
  cat(
    "\nIndividual plots:\n",
    cosinor_dir,
    "\n"
  )
  
}


############################
# 35. Save DESeq2 object
############################

saveRDS(
  
  dds,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_dds.rds"
  )
)


############################
# 36. Final summary
############################

cat(
  "\n\n============================================\n"
)

cat(
  "DESeq2 24-hour rhythmic analysis COMPLETE\n"
)

cat(
  "============================================\n\n"
)


cat(
  "Genes tested after filtering: ",
  nrow(dds),
  "\n",
  sep = ""
)


cat(
  "FDR < 0.05 rhythmic genes: ",
  nrow(rhythmic_results),
  "\n",
  sep = ""
)


cat(
  "Top genes plotted: ",
  length(top_genes),
  "\n",
  sep = ""
)


cat(
  "\nOutput directory:\n",
  outdir,
  "\n"
)


############################
# 37. Print top genes
############################

if (nrow(rhythmic_results) > 0) {
  
  cat(
    "\nTop rhythmic genes by FDR:\n\n"
  )
  
  
  columns_to_print <- c(
    
    "baseMean",
    
    "pvalue",
    
    "padj",
    
    "beta_cos",
    
    "beta_sin",
    
    "amplitude",
    
    "peak_ZT"
  )
  
  
  # Only use columns that actually exist
  
  columns_to_print <- columns_to_print[
    columns_to_print %in%
      colnames(rhythmic_results)
  ]
  
  
  print(
    
    rhythmic_results[
      
      1:min(
        20,
        nrow(rhythmic_results)
      ),
      
      columns_to_print,
      
      drop = FALSE
    ]
  )
  
  
} else {
  
  cat(
    "\nNo genes passed FDR < 0.05.\n"
  )
  
}


############################
# 38. Final message
############################

cat(
  "\nAnalysis finished successfully.\n"
)

cat(
  "Results written to:\n",
  outdir,
  "\n"
)

############################################################
# END OF SCRIPT
############################################################
