# ============================================================
# Target clock gene plots from existing DESeq2 objects
# ============================================================

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# Check required objects from the main analysis
# ------------------------------------------------------------

required_objects <- c(
  "count_matrix",
  "outdir",
  "sample_info",
  "dds",
  "norm_counts",
  "beta_cos",
  "beta_sin",
  "amplitude",
  "phase_hours",
  "omega",
  "rhythm_results"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing from the main analysis:\n",
    paste(missing_objects, collapse = ", "),
    "\n\nRun the main DESeq2 script first."
  )
}

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

target_outdir <- file.path(outdir, "target_clock_genes")

if (!dir.exists(target_outdir)) {
  dir.create(target_outdir, recursive = TRUE)
}

# ------------------------------------------------------------
# Target clock genes and aliases
# ------------------------------------------------------------

target_genes <- c(
  "CLOCK",
  "BMAL1",
  "PER1",
  "PER2",
  "PER3",
  "CRY1",
  "CRY2",
  "REV-ERB",
  "ROR",
  "DBP",
  "NFIL3"
)

target_aliases <- list(
  "CLOCK"   = c("CLOCK"),
  "BMAL1"   = c("BMAL1", "ARNTL"),
  "PER1"    = c("PER1"),
  "PER2"    = c("PER2"),
  "PER3"    = c("PER3"),
  "CRY1"    = c("CRY1"),
  "CRY2"    = c("CRY2"),
  "REV-ERB" = c("REV-ERB", "REV-ERBA", "REV-ERBB", "NR1D1", "NR1D2"),
  "ROR"     = c("ROR", "RORA", "RORB", "RORC"),
  "DBP"     = c("DBP"),
  "NFIL3"   = c("NFIL3")
)

# ------------------------------------------------------------
# Validate sample information
# ------------------------------------------------------------

if (!"ZT" %in% colnames(sample_info)) {
  stop("sample_info must contain a column named 'ZT'.")
}

if (!"replicate" %in% colnames(sample_info)) {
  stop("sample_info must contain a column named 'replicate'.")
}

if (is.null(rownames(sample_info))) {
  stop("sample_info must have sample names as row names.")
}

if (!all(rownames(sample_info) %in% colnames(count_matrix))) {
  stop("Not all sample_info row names are present in count_matrix.")
}

if (!all(rownames(sample_info) %in% colnames(norm_counts))) {
  stop("Not all sample_info row names are present in norm_counts.")
}

# ------------------------------------------------------------
# Match target genes case-insensitively
# ------------------------------------------------------------

available_gene_ids <- trimws(rownames(count_matrix))
available_gene_ids_lower <- tolower(available_gene_ids)

target_matches <- list()

for (target in target_genes) {
  
  aliases <- target_aliases[[target]]
  aliases_lower <- tolower(trimws(aliases))
  
  matches <- available_gene_ids[
    available_gene_ids_lower %in% aliases_lower
  ]
  
  target_matches[[target]] <- matches
}

# ------------------------------------------------------------
# Print target gene matches
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Target clock gene lookup\n")
cat("============================================================\n\n")

for (target in target_genes) {
  
  matches <- target_matches[[target]]
  
  if (length(matches) == 0) {
    
    cat(target, ": NOT FOUND\n")
    
  } else {
    
    cat(
      target,
      ": ",
      paste(matches, collapse = ", "),
      "\n",
      sep = ""
    )
  }
}

cat("\n")

# ------------------------------------------------------------
# Stop if no target genes were found
# ------------------------------------------------------------

if (!any(lengths(target_matches) > 0)) {
  
  cat("The first Geneid values are:\n")
  print(head(available_gene_ids, 20))
  
  stop(
    "None of the target genes were found in the Geneid values of the featureCounts table."
  )
}

# ------------------------------------------------------------
# Save target gene lookup table
# ------------------------------------------------------------

target_lookup <- data.frame(
  Target = target_genes,
  Geneid = vapply(
    target_genes,
    function(x) {
      matches <- target_matches[[x]]
      if (length(matches) == 0) NA_character_ else matches[1]
    },
    character(1)
  ),
  All_matches = vapply(
    target_genes,
    function(x) {
      matches <- target_matches[[x]]
      if (length(matches) == 0) {
        NA_character_
      } else {
        paste(matches, collapse = "; ")
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  target_lookup,
  file.path(target_outdir, "target_gene_matches.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Select first matching Geneid for each target
# ------------------------------------------------------------

selected_targets <- target_lookup[
  !is.na(target_lookup$Geneid),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Save raw target counts
# ------------------------------------------------------------

selected_gene_ids <- selected_targets$Geneid

raw_target_counts <- count_matrix[
  selected_gene_ids,
  ,
  drop = FALSE
]

rownames(raw_target_counts) <- selected_targets$Target

write.csv(
  raw_target_counts,
  file.path(target_outdir, "target_clock_gene_raw_counts.csv")
)

# ------------------------------------------------------------
# Save raw target counts with original Geneid
# ------------------------------------------------------------

raw_target_counts_labelled <- data.frame(
  Target = selected_targets$Target,
  Geneid = selected_targets$Geneid,
  raw_target_counts,
  check.names = FALSE
)

write.csv(
  raw_target_counts_labelled,
  file.path(
    target_outdir,
    "target_clock_gene_raw_counts_labelled.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Long-format raw target counts
# ------------------------------------------------------------

raw_target_long <- do.call(
  rbind,
  lapply(seq_len(nrow(selected_targets)), function(i) {
    
    gene_id <- selected_targets$Geneid[i]
    target_label <- selected_targets$Target[i]
    
    data.frame(
      Target = target_label,
      Geneid = gene_id,
      sample = colnames(count_matrix),
      raw_count = as.numeric(count_matrix[gene_id, ]),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  raw_target_long,
  file.path(
    target_outdir,
    "target_clock_gene_raw_counts_long.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# DESeq2 coefficient matrix
# ------------------------------------------------------------

coef_matrix <- coef(dds)

intercept_name <- grep(
  "Intercept",
  colnames(coef_matrix),
  ignore.case = TRUE,
  value = TRUE
)[1]

if (is.na(intercept_name) || length(intercept_name) == 0) {
  stop("Could not identify the intercept coefficient in coef(dds).")
}

if (!"cos24" %in% colnames(coef_matrix)) {
  stop("Could not find 'cos24' in coef(dds).")
}

if (!"sin24" %in% colnames(coef_matrix)) {
  stop("Could not find 'sin24' in coef(dds).")
}

# ------------------------------------------------------------
# Function to retrieve FDR / adjusted p-value
# ------------------------------------------------------------

get_padj <- function(gene_id) {
  
  if ("gene_id" %in% colnames(rhythm_results)) {
    
    idx <- which(
      rhythm_results$gene_id == gene_id
    )
    
    if (length(idx) > 0) {
      return(rhythm_results$padj[idx[1]])
    }
  }
  
  if (!is.null(rownames(rhythm_results))) {
    
    idx <- which(
      rownames(rhythm_results) == gene_id
    )
    
    if (length(idx) > 0) {
      return(rhythm_results$padj[idx[1]])
    }
  }
  
  return(NA_real_)
}

# ------------------------------------------------------------
# Plot function
# ------------------------------------------------------------

plot_target_gene <- function(target_label, gene_id) {
  
  gene_counts <- norm_counts[
    gene_id,
    ,
    drop = TRUE
  ]
  
  plot_data <- data.frame(
    sample = names(gene_counts),
    ZT = sample_info[names(gene_counts), "ZT"],
    replicate = sample_info[
      names(gene_counts),
      "replicate"
    ],
    expression = as.numeric(gene_counts),
    stringsAsFactors = FALSE
  )
  
  mean_data <- aggregate(
    expression ~ ZT,
    data = plot_data,
    FUN = mean
  )
  
  sd_data <- aggregate(
    expression ~ ZT,
    data = plot_data,
    FUN = sd
  )
  
  summary_data <- merge(
    mean_data,
    sd_data,
    by = "ZT",
    suffixes = c("", "_SD")
  )
  
  colnames(summary_data) <- c(
    "ZT",
    "expression",
    "SD"
  )
  
  # ----------------------------------------------------------
  # Harmonic coefficients
  # ----------------------------------------------------------
  
  gene_intercept <- coef_matrix[
    gene_id,
    intercept_name
  ]
  
  gene_beta_cos <- coef_matrix[
    gene_id,
    "cos24"
  ]
  
  gene_beta_sin <- coef_matrix[
    gene_id,
    "sin24"
  ]
  
  # ----------------------------------------------------------
  # Smooth fitted curve
  # ----------------------------------------------------------
  
  prediction_ZT <- seq(
    0,
    24,
    length.out = 500
  )
  
  prediction_log2 <-
    gene_intercept +
    gene_beta_cos * cos(omega * prediction_ZT) +
    gene_beta_sin * sin(omega * prediction_ZT)
  
  prediction_fit <- 2^prediction_log2
  
  prediction_data <- data.frame(
    ZT = prediction_ZT,
    fit = prediction_fit
  )
  
  # ----------------------------------------------------------
  # Rhythm statistics
  # ----------------------------------------------------------
  
  gene_peak <- phase_hours[gene_id]
  
  gene_amplitude <- amplitude[gene_id]
  
  gene_padj <- get_padj(gene_id)
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    plot_data,
    aes(
      x = ZT,
      y = expression
    )
  ) +
    
    geom_jitter(
      width = 0.35,
      height = 0,
      size = 2.4,
      alpha = 0.65
    ) +
    
    geom_errorbar(
      data = summary_data,
      aes(
        x = ZT,
        ymin = pmax(0, expression - SD),
        ymax = expression + SD
      ),
      width = 0.4,
      linewidth = 0.7
    ) +
    
    geom_point(
      data = summary_data,
      aes(
        x = ZT,
        y = expression
      ),
      size = 4
    ) +
    
    geom_line(
      data = prediction_data,
      aes(
        x = ZT,
        y = fit
      ),
      linewidth = 1.2
    ) +
    
    geom_vline(
      xintercept = gene_peak,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
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
      limits = c(0, 24)
    ) +
    
    expand_limits(
      y = 0
    ) +
    
    labs(
      title = target_label,
      subtitle = paste0(
        "FDR = ",
        format.pval(
          gene_padj,
          digits = 3,
          eps = 1e-10
        ),
        "   |   Amplitude = ",
        round(gene_amplitude, 3),
        " log2 units",
        "   |   Peak ZT = ",
        round(gene_peak, 2)
      ),
      x = "Zeitgeber Time (ZT)",
      y = "DESeq2-normalized expression"
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      axis.title = element_text(
        face = "bold"
      )
    )
  
  return(p)
}

# ------------------------------------------------------------
# Generate individual target plots
# ------------------------------------------------------------

target_plots <- list()

for (i in seq_len(nrow(selected_targets))) {
  
  target_label <- selected_targets$Target[i]
  gene_id <- selected_targets$Geneid[i]
  
  cat(
    "Generating plot for ",
    target_label,
    " (",
    gene_id,
    ")\n",
    sep = ""
  )
  
  p <- plot_target_gene(
    target_label = target_label,
    gene_id = gene_id
  )
  
  target_plots[[target_label]] <- p
  
  ggsave(
    file.path(
      target_outdir,
      paste0(
        target_label,
        "_rhythm.pdf"
      )
    ),
    p,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(
      target_outdir,
      paste0(
        target_label,
        "_rhythm.png"
      )
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ------------------------------------------------------------
# Figure 1: Core clock genes
# ------------------------------------------------------------

figure1_available <- intersect(
  c(
    "CLOCK",
    "BMAL1",
    "PER1",
    "PER2",
    "PER3",
    "CRY1"
  ),
  names(target_plots)
)

if (length(figure1_available) > 0) {
  
  p <- wrap_plots(
    target_plots[figure1_available],
    ncol = 2,
    nrow = 3
  )
  
  ggsave(
    file.path(
      target_outdir,
      "Figure1_Core_Clock_Genes.pdf"
    ),
    p,
    width = 12,
    height = 16
  )
  
  ggsave(
    file.path(
      target_outdir,
      "Figure1_Core_Clock_Genes.png"
    ),
    p,
    width = 12,
    height = 16,
    dpi = 300
  )
}

# ------------------------------------------------------------
# Figure 2: Additional clock genes
# ------------------------------------------------------------

figure2_available <- intersect(
  c(
    "CRY2",
    "REV-ERB",
    "ROR",
    "DBP",
    "NFIL3"
  ),
  names(target_plots)
)

if (length(figure2_available) > 0) {
  
  p <- wrap_plots(
    target_plots[figure2_available],
    ncol = 2,
    nrow = 3
  )
  
  ggsave(
    file.path(
      target_outdir,
      "Figure2_Additional_Clock_Genes.pdf"
    ),
    p,
    width = 12,
    height = 16
  )
  
  ggsave(
    file.path(
      target_outdir,
      "Figure2_Additional_Clock_Genes.png"
    ),
    p,
    width = 12,
    height = 16,
    dpi = 300
  )
}

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

target_summary <- do.call(
  rbind,
  lapply(
    seq_len(nrow(selected_targets)),
    function(i) {
      
      target_label <- selected_targets$Target[i]
      gene_id <- selected_targets$Geneid[i]
      
      data.frame(
        Target = target_label,
        Geneid = gene_id,
        Intercept = coef_matrix[
          gene_id,
          intercept_name
        ],
        beta_cos = beta_cos[gene_id],
        beta_sin = beta_sin[gene_id],
        amplitude = amplitude[gene_id],
        Peak_ZT = phase_hours[gene_id],
        FDR = get_padj(gene_id),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  target_summary,
  file.path(
    target_outdir,
    "target_clock_gene_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Save RDS objects
# ------------------------------------------------------------

saveRDS(
  target_matches,
  file.path(
    target_outdir,
    "target_matches.rds"
  )
)

saveRDS(
  target_plots,
  file.path(
    target_outdir,
    "target_plots.rds"
  )
)

saveRDS(
  target_summary,
  file.path(
    target_outdir,
    "target_summary.rds"
  )
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Target clock gene analysis complete\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  target_outdir,
  "\n\n"
)

cat(
  "Target genes found:\n",
  paste(
    selected_targets$Target,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Summary table:\n",
  file.path(
    target_outdir,
    "target_clock_gene_summary.csv"
  ),
  "\n\n"
)

cat(
  "Figure 1:\n",
  file.path(
    target_outdir,
    "Figure1_Core_Clock_Genes.pdf"
  ),
  "\n\n"
)

cat(
  "Figure 2:\n",
  file.path(
    target_outdir,
    "Figure2_Additional_Clock_Genes.pdf"
  ),
  "\n\n"
)
