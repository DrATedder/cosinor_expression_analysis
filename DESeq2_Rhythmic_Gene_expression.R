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

library(DESeq2)
library(ggplot2)
library(pheatmap)


############################
# 2. Input/output settings
############################

count_file <- "/home/andrew/NMR_featureCounts_clean.txt"

outdir <- "/home/andrew/DESeq2_rhythm"

# Number of significant genes to plot
TOP_N <- 50


dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################
# 3. Read featureCounts data
############################

cat("\nReading featureCounts data...\n")

counts_raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
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
# 4. Extract count matrix
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

count_matrix <- as.matrix(
  count_matrix
)

storage.mode(count_matrix) <- "integer"


############################
# 5. Sample metadata
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
# 6. Check sample names
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


cat("Sample names match successfully.\n")


############################
# 7. Create harmonic variables
############################

omega <- 2 * pi / 24


sample_info$cos24 <- cos(
  omega * sample_info$ZT
)


sample_info$sin24 <- sin(
  omega * sample_info$ZT
)


############################
# 8. Save metadata
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
# 9. Filter low-count genes
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


count_matrix_filtered <- count_matrix[
  keep,
  ,
  drop = FALSE
]


############################
# 10. Create DESeqDataSet
############################

dds <- DESeqDataSetFromMatrix(
  
  countData = count_matrix_filtered,
  
  colData = sample_info,
  
  design = ~ cos24 + sin24
)


############################
# 11. Run DESeq2 LRT
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
# 12. Display model
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
# 13. Extract LRT results
############################

res <- results(
  dds
)


res <- res[
  order(res$padj),
]


############################
# 14. Save ALL results
############################

write.csv(
  
  as.data.frame(res),
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_all_genes.csv"
  )
)


############################
# 15. Select significant genes
############################

rhythmic <- res[
  
  !is.na(res$padj) &
    
    res$padj < 0.05,
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
# 16. Save FDR < 0.05 results
############################

write.csv(
  
  as.data.frame(rhythmic),
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythmic_genes_FDR05.csv"
  )
)


############################
# 17. Extract DESeq2 coefficients
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


if (!"cos24" %in% colnames(coef_matrix)) {
  
  stop(
    "cos24 coefficient not found."
  )
}


if (!"sin24" %in% colnames(coef_matrix)) {
  
  stop(
    "sin24 coefficient not found."
  )
}


beta_cos <- coef_matrix[
  ,
  "cos24"
]


beta_sin <- coef_matrix[
  ,
  "sin24"
]


############################
# 18. Calculate amplitude
############################

amplitude <- sqrt(
  
  beta_cos^2 +
    
    beta_sin^2
)


############################
# 19. Calculate peak phase
############################

phase_radians <- atan2(
  
  beta_sin,
  
  beta_cos
)


phase_hours <- (
  
  phase_radians %% (2 * pi)
  
) * 24 / (2 * pi)


############################
# 20. Create complete rhythm table
############################

rhythm_results <- as.data.frame(
  res
)


rhythm_results$beta_cos <- beta_cos[
  match(
    rownames(rhythm_results),
    names(beta_cos)
  )
]


rhythm_results$beta_sin <- beta_sin[
  match(
    rownames(rhythm_results),
    names(beta_sin)
  )
]


rhythm_results$amplitude <- amplitude[
  match(
    rownames(rhythm_results),
    names(amplitude)
  )
]


rhythm_results$peak_ZT <- phase_hours[
  match(
    rownames(rhythm_results),
    names(phase_hours)
  )
]


rhythm_results <- rhythm_results[
  order(rhythm_results$padj),
]


############################
# 21. Save complete rhythm table
############################

write.csv(
  
  rhythm_results,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_with_phase_amplitude.csv"
  )
)


############################
# 22. Significant rhythm table
############################

rhythmic_results <- rhythm_results[
  
  !is.na(rhythm_results$padj) &
    
    rhythmic_results$padj < 0.05,
  
]


rhythmic_results <- rhythmic_results[
  order(rhythmic_results$padj),
]


write.csv(
  
  rhythmic_results,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv"
  )
)


############################
# 23. Normalized counts
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
  )
)


############################
# 24. VST
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
  )
)


############################
# 25. PCA
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
# 26. Sample correlation heatmap
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
    sample_info[
      ,
      "ZT",
      drop = FALSE
    ],
  
  annotation_row =
    sample_info[
      ,
      "ZT",
      drop = FALSE
    ],
  
  main = "Sample correlation"
)


dev.off()


############################
# 27. MA plot
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
# 28. Top 20 rhythmic gene heatmap
############################

top_heatmap_genes <- rownames(
  
  rhythmic_results[
    1:min(
      20,
      nrow(rhythmic_results)
    ),
    ,
    drop = FALSE
  ]
)


if (length(top_heatmap_genes) > 0) {
  
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
}


############################################################
# 29. COSINOR PLOTS
#
# IMPORTANT:
#
# Only the TOP_N genes from the FDR < 0.05 table
# are plotted.
#
# They are ranked by lowest adjusted p-value.
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
  
} else {
  
  
  ##########################################################
  # 29.1 Select top genes
  ##########################################################
  
  top_gene_count <- min(
    
    TOP_N,
    
    nrow(rhythmic_results)
  )
  
  
  top_genes <- rownames(
    
    rhythmic_results[
      1:top_gene_count,
      ,
      drop = FALSE
    ]
  )
  
  
  top_results <- rhythmic_results[
    
    top_genes,
    
    ,
    drop = FALSE
  ]
  
  
  ##########################################################
  # 29.2 Save top-gene table
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
    )
  )
  
  
  ##########################################################
  # 29.3 Cosinor output directory
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
  # 29.4 Combined PDF
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
  # 29.5 Check fitted means
  ##########################################################
  
  if (!"mu" %in% assayNames(dds)) {
    
    dev.off()
    
    stop(
      "DESeq2 fitted means ('mu') are not available."
    )
  }
  
  
  fitted_mu_matrix <- assays(dds)[[
    "mu"
  ]]
  
  
  ##########################################################
  # 29.6 Loop through top genes
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
    # SD
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
    
    
    ########################################################
    # DESeq2 fitted means
    ########################################################
    
    fitted_mu <- fitted_mu_matrix[
      
      gene,
      
      ,
      
      drop = TRUE
    ]
    
    
    ########################################################
    # Size factors
    ########################################################
    
    gene_size_factors <- sizeFactors(
      dds
    )
    
    
    ########################################################
    # Convert fitted counts to normalized expression
    ########################################################
    
    fitted_q <- fitted_mu /
      gene_size_factors
    
    
    ########################################################
    # Recover intercept
    #
    # This avoids assuming that the coefficient is called
    # "(Intercept)".
    ########################################################
    
    intercept_estimates <- (
      
      log2(
        pmax(
          fitted_q,
          1e-8
        )
      )
      
      -
        
        gene_beta_cos *
        sample_info$cos24
      
      -
        
        gene_beta_sin *
        sample_info$sin24
    )
    
    
    finite_intercepts <-
      intercept_estimates[
        is.finite(
          intercept_estimates
        )
      ]
    
    
    if (length(finite_intercepts) == 0) {
      
      warning(
        "Could not calculate intercept for ",
        gene,
        ". Skipping."
      )
      
      next
    }
    
    
    gene_intercept <- mean(
      finite_intercepts
    )
    
    
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
    # Convert back to normalized expression scale
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
      gene,
      "padj"
    ]
    
    
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
# 30. Save DESeq2 object
############################

saveRDS(
  
  dds,
  
  file = file.path(
    outdir,
    "DESeq2_24h_rhythm_dds.rds"
  )
)


############################
# 31. Final summary
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


if (exists("top_genes")) {
  
  cat(
    "Top genes plotted: ",
    length(top_genes),
    "\n",
    sep = ""
  )
}


cat(
  "\nOutput directory:\n",
  outdir,
  "\n"
)


if (nrow(rhythmic_results) > 0) {
  
  cat(
    "\nTop rhythmic genes by FDR:\n\n"
  )
  
  
  print(
    
    rhythmic_results[
      
      1:min(
        20,
        nrow(rhythmic_results)
      ),
      
      c(
        "baseMean",
        "pvalue",
        "padj",
        "beta_cos",
        "beta_sin",
        "amplitude",
        "peak_ZT"
      ),
      
      drop = FALSE
    ]
  )
}


cat(
  "\nAnalysis finished successfully.\n"
)
