############################################################
# COSINOR PLOTS COLOURED BY GONADAL SEX
#
# Uses outputs from the original DESeq2 analysis.
#
#
# Individual replicate points:
#   Coloured by gonadal sex (M/F)
#
#
# Cosinor curve:
#   Original pooled DESeq2 model fit.
#
# Outputs:
#   top_50_cosinor_by_sex.pdf
#   top_50_cosinor_by_sex/
#       Individual PDF and PNG plots
############################################################


############################
# 1. Load packages
############################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})


############################
# 2. Input/output settings
############################

outdir <- "/home/andrew/DESeq2_rhythm"

TOP_N <- 50

cosinor_dir <- file.path(
  outdir,
  "top_50_cosinor_by_sex"
)

dir.create(
  cosinor_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################
# 3. Load existing results
############################

# Normalized counts from the original analysis

norm_counts <- read.csv(
  file.path(
    outdir,
    "DESeq2_normalized_counts.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

# Rhythm results, including FDR, amplitude and peak ZT

rhythm_results <- read.csv(
  file.path(
    outdir,
    "DESeq2_24h_rhythm_with_phase_amplitude.csv"
  ),
  check.names = FALSE
)

# Original fitted DESeq2 model

dds <- readRDS(
  file.path(
    outdir,
    "DESeq2_24h_rhythm_dds.rds"
  )
)

coef_matrix <- coef(dds)

coef_names <- colnames(coef_matrix)


############################
# 4. Identify intercept
############################

intercept_name <- coef_names[
  grepl(
    "Intercept",
    coef_names,
    ignore.case = TRUE
  )
]

if (length(intercept_name) != 1) {
  stop(
    "Could not uniquely identify the DESeq2 intercept."
  )
}

intercept_name <- intercept_name[1]


############################
# 5. Sample metadata
############################

# Preserve the original sample-to-ZT assignments.

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
  
  ZT = rep(
    c(0, 4, 8, 12, 16, 20),
    each = 4
  ),
  
  replicate = rep(
    1:4,
    times = 6
  ),
  
  stringsAsFactors = FALSE
)


############################
# 6. Gonadal sex information
############################

sex_info <- data.frame(
  
  sample = sprintf(
    "NMR2-COL1-HYP-%02d",
    1:24
  ),
  
  sex = c(
    "M", "M", "M", "F",
    "M", "M", "M", "M",
    "F", "F", "M", "F",
    "F", "M", "M", "F",
    "F", "M", "F", "F",
    "M", "F", "F", "F"
  ),
  
  stringsAsFactors = FALSE
)


############################
# 7. Merge metadata
############################

sample_info$sex <- sex_info$sex[
  match(
    sample_info$sample,
    sex_info$sample
  )
]

if (anyNA(sample_info$sex)) {
  stop(
    "Missing gonadal sex assignments detected."
  )
}

sample_info$sex <- factor(
  sample_info$sex,
  levels = c("F", "M")
)

rownames(sample_info) <- sample_info$sample


############################
# 8. Validate input data
############################

if (!all(
  sample_info$sample %in% colnames(norm_counts)
)) {
  stop(
    "Some metadata samples are missing from normalized counts."
  )
}

if (!all(
  sample_info$sample %in% colnames(dds)
)) {
  stop(
    "Some metadata samples are missing from the DESeq2 object."
  )
}

required_columns <- c(
  "gene_id",
  "padj",
  "amplitude",
  "peak_ZT"
)

if (!all(
  required_columns %in% colnames(rhythm_results)
)) {
  stop(
    "Required columns are missing from the rhythm results file."
  )
}

if (!all(
  c("cos24", "sin24") %in% coef_names
)) {
  stop(
    "The cos24 or sin24 coefficient is missing."
  )
}


############################
# 9. Select top 50 rhythmic genes
############################

rhythmic_results <- rhythm_results[
  !is.na(rhythm_results$padj) &
    rhythm_results$padj < 0.05,
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

top_gene_count <- min(
  TOP_N,
  nrow(rhythmic_results)
)

if (top_gene_count == 0) {
  stop(
    "No genes passed FDR < 0.05. No plots can be generated."
  )
}

top_genes <- head(
  rhythmic_results$gene_id,
  top_gene_count
)

cat(
  "\nGenerating cosinor plots for ",
  length(top_genes),
  " genes.\n",
  sep = ""
)


############################
# 10. Prepare cosinor predictions
############################

omega <- 2 * pi / 24

prediction_ZT <- seq(
  0,
  24,
  length.out = 500
)


############################
# 11. Create combined PDF
############################

combined_pdf <- file.path(
  outdir,
  "top_50_cosinor_by_sex.pdf"
)

pdf(
  combined_pdf,
  width = 8,
  height = 6
)


############################
# 12. Generate individual plots
############################

for (gene in top_genes) {
  
  cat(
    "Plotting: ",
    gene,
    "\n",
    sep = ""
  )
  
  
  ############################
  # 12.1 Observed expression
  ############################
  
  gene_counts <- norm_counts[
    gene,
    sample_info$sample,
    drop = TRUE
  ]
  
  plot_data <- data.frame(
    
    sample = sample_info$sample,
    
    ZT = sample_info$ZT,
    
    replicate = sample_info$replicate,
    
    sex = sample_info$sex,
    
    expression = as.numeric(
      gene_counts
    )
  )
  
  
  ############################
  # 12.2 Original mean expression
  ############################
  
  # Mean calculated across all four replicates,
  # irrespective of sex.
  
  mean_data <- aggregate(
    expression ~ ZT,
    data = plot_data,
    FUN = mean
  )
  
  
  ############################
  # 12.3 Original standard deviation
  ############################
  
  # SD calculated across all four replicates,
  # irrespective of sex.
  
  sd_data <- aggregate(
    expression ~ ZT,
    data = plot_data,
    FUN = sd
  )
  
  names(sd_data)[2] <- "SD"
  
  
  ############################
  # 12.4 Original summary data
  ############################
  
  summary_data <- merge(
    mean_data,
    sd_data,
    by = "ZT"
  )
  
  
  ############################
  # 12.5 Retrieve DESeq2 coefficients
  ############################
  
  gene_beta_cos <- coef_matrix[
    gene,
    "cos24"
  ]
  
  gene_beta_sin <- coef_matrix[
    gene,
    "sin24"
  ]
  
  gene_intercept <- coef_matrix[
    gene,
    intercept_name
  ]
  
  
  ############################
  # 12.6 Check coefficients
  ############################
  
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
  
  
  ############################
  # 12.7 Original cosinor fit
  ############################
  
  prediction_log2 <- (
    
    gene_intercept +
      
      gene_beta_cos *
      cos(omega * prediction_ZT) +
      
      gene_beta_sin *
      sin(omega * prediction_ZT)
    
  )
  
  prediction_fit <- 2^prediction_log2
  
  prediction_data <- data.frame(
    
    ZT = prediction_ZT,
    
    fit = prediction_fit
  )
  
  
  ############################
  # 12.8 Rhythm statistics
  ############################
  
  gene_info <- rhythmic_results[
    rhythmic_results$gene_id == gene,
    ,
    drop = FALSE
  ]
  
  gene_padj <- gene_info$padj[1]
  
  gene_amplitude <- gene_info$amplitude[1]
  
  gene_peak <- gene_info$peak_ZT[1]
  
  
  ############################
  # 12.9 Plot
  ############################
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = ZT,
      y = expression
    )
    
  ) +
    
    
    ############################
  # Individual replicate points
  ############################
  
  # Sex colouring is applied ONLY to these points.
  
  geom_jitter(
    
    aes(
      colour = sex
    ),
    
    width = 0.35,
    
    height = 0,
    
    size = 2.4,
    
    alpha = 0.65
  ) +
    
    
    ############################
  # Original error bars
  ############################
  
  # Unchanged: mean +/- SD across all four replicates.
  
  geom_errorbar(
    
    data = summary_data,
    
    aes(
      
      x = ZT,
      
      ymin = pmax(
        0,
        expression - SD
      ),
      
      ymax = expression + SD
    ),
    
    width = 0.4,
    
    linewidth = 0.7,
    
    inherit.aes = FALSE
  ) +
    
    
    ############################
  # Original mean points
  ############################
  
  # Unchanged: one mean point per ZT.
  
  geom_point(
    
    data = summary_data,
    
    aes(
      
      x = ZT,
      
      y = expression
    ),
    
    size = 4,
    
    inherit.aes = FALSE
  ) +
    
    
    ############################
  # Original cosinor fitted curve
  ############################
  
  geom_line(
    
    data = prediction_data,
    
    aes(
      
      x = ZT,
      
      y = fit
    ),
    
    linewidth = 1.2,
    
    inherit.aes = FALSE
  ) +
    
    
    ############################
  # Original predicted peak
  ############################
  
  geom_vline(
    
    xintercept = gene_peak,
    
    linetype = "dashed",
    
    linewidth = 0.7
  ) +
    
    
    ############################
  # Sex colours for individual points
  ############################
  
  scale_colour_manual(
    
    values = c(
      "F" = "#D55E00",
      "M" = "#0072B2"
    ),
    
    labels = c(
      "F" = "Female",
      "M" = "Male"
    ),
    
    name = "Gonadal sex"
  ) +
    
    
    ############################
  # Original X axis
  ############################
  
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
    
    
    ############################
  # Original Y axis
  ############################
  
  expand_limits(
    y = 0
  ) +
    
    
    ############################
  # Original labels
  ############################
  
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
    
    
    ############################
  # Original theme
  ############################
  
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
  
  
  ############################
  # 12.10 Add to combined PDF
  ############################
  
  print(p)
  
  
  ############################
  # 12.11 Safe filename
  ############################
  
  safe_gene_name <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    gene
  )
  
  
  ############################
  # 12.12 Individual PDF
  ############################
  
  ggsave(
    
    filename = file.path(
      
      cosinor_dir,
      
      paste0(
        safe_gene_name,
        "_cosinor_by_sex.pdf"
      )
    ),
    
    plot = p,
    
    width = 8,
    
    height = 6
  )
  
  
  ############################
  # 12.13 Individual PNG
  ############################
  
  ggsave(
    
    filename = file.path(
      
      cosinor_dir,
      
      paste0(
        safe_gene_name,
        "_cosinor_by_sex.png"
      )
    ),
    
    plot = p,
    
    width = 8,
    
    height = 6,
    
    dpi = 300
  )
  
}


############################
# 13. Close combined PDF
############################

dev.off()


############################
# 14. Final summary
############################

cat(
  "\n============================================\n"
)

cat(
  "SEX-COLOURED COSINOR PLOTTING COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "Significant rhythmic genes available: ",
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
  "\nCombined PDF:\n",
  combined_pdf,
  "\n"
)

cat(
  "\nIndividual plots:\n",
  cosinor_dir,
  "\n"
)

############################################################
# END
############################################################
