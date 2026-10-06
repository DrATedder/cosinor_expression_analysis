############################################################
# Biological Interpretation of 24-h Rhythmic Genes
# Species: Heterocephalus glaber (naked mole-rat)
#
# Input:
#   DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv
#
# Includes:
#   - Rhythmic gene summary
#   - Phase distribution
#   - Phase bins
#   - Amplitude distribution
#   - Phase vs amplitude
#   - Phase-ordered VST heatmap
#   - g:Profiler GO Biological Process enrichment
#   - g:Profiler KEGG enrichment
#   - Overall enrichment
#   - Phase-specific enrichment
#   - Diagnostics for recognised/unrecognised gene IDs
#
# IMPORTANT:
#   Gene IDs are taken explicitly from the FIRST column.
#   No automatic gene-ID column detection is used.
############################################################


############################
# 1. Install/load packages
############################

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "pheatmap",
  "RColorBrewer",
  "gprofiler2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(gprofiler2)


############################
# 2. File paths
############################

input_file <- "/home/andrew/DESeq2_rhythm/DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv"

output_dir <- "/home/andrew/DESeq2_rhythm/biological_interpretation"

vst_file <- "/home/andrew/DESeq2_rhythm/DESeq2_VST_expression.csv"

background_file <- "/home/andrew/DESeq2_rhythm/DESeq2_24h_rhythm_all_genes.csv"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################
# 3. Analysis settings
############################

# Species:
# Heterocephalus glaber = naked mole-rat
gprofiler_organism <- "hgfemale"

# g:Profiler sources available/relevant here
gprofiler_sources <- c(
  "GO:BP",
  "KEGG"
)

# Rhythmicity threshold
padj_threshold <- 0.05

# Minimum number of genes required for
# phase-specific enrichment
minimum_phase_genes <- 5


############################
# 4. Read rhythmic results
############################

cat("\n============================================\n")
cat("READING RHYTHMIC GENE RESULTS\n")
cat("============================================\n\n")

rhythmic <- read_csv(
  input_file,
  show_col_types = FALSE
)

# Explicitly define first column as gene ID
colnames(rhythmic)[1] <- "gene"

rhythmic$gene <- as.character(rhythmic$gene)

cat("Number of rows in input file:", nrow(rhythmic), "\n")
cat("Number of columns:", ncol(rhythmic), "\n\n")

cat("Columns detected:\n")
print(colnames(rhythmic))


############################
# 5. Basic filtering
############################

rhythmic_filtered <- rhythmic %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(padj),
    padj < padj_threshold,
    !is.na(peak_ZT),
    !is.na(amplitude)
  )

cat("\n============================================\n")
cat("RHYTHMIC GENE FILTERING\n")
cat("============================================\n\n")

cat(
  "Rhythmic genes with padj <",
  padj_threshold,
  ":",
  nrow(rhythmic_filtered),
  "\n"
)


############################
# 6. Create phase bins
############################

rhythmic_filtered <- rhythmic_filtered %>%
  mutate(
    peak_ZT = as.numeric(peak_ZT),
    amplitude = as.numeric(amplitude),
    phase_bin = case_when(
      peak_ZT >= 0  & peak_ZT < 4  ~ "ZT0-4",
      peak_ZT >= 4  & peak_ZT < 8  ~ "ZT4-8",
      peak_ZT >= 8  & peak_ZT < 12 ~ "ZT8-12",
      peak_ZT >= 12 & peak_ZT < 16 ~ "ZT12-16",
      peak_ZT >= 16 & peak_ZT < 20 ~ "ZT16-20",
      peak_ZT >= 20 & peak_ZT <= 24 ~ "ZT20-24",
      TRUE ~ NA_character_
    )
  )

phase_levels <- c(
  "ZT0-4",
  "ZT4-8",
  "ZT8-12",
  "ZT12-16",
  "ZT16-20",
  "ZT20-24"
)

rhythmic_filtered$phase_bin <- factor(
  rhythmic_filtered$phase_bin,
  levels = phase_levels
)


############################
# 7. Summary statistics
############################

summary_table <- data.frame(
  statistic = c(
    "Number of rhythmic genes",
    "Minimum adjusted p-value",
    "Median adjusted p-value",
    "Maximum adjusted p-value",
    "Minimum amplitude",
    "Median amplitude",
    "Maximum amplitude",
    "Minimum peak ZT",
    "Median peak ZT",
    "Maximum peak ZT"
  ),
  value = c(
    nrow(rhythmic_filtered),
    min(rhythmic_filtered$padj, na.rm = TRUE),
    median(rhythmic_filtered$padj, na.rm = TRUE),
    max(rhythmic_filtered$padj, na.rm = TRUE),
    min(rhythmic_filtered$amplitude, na.rm = TRUE),
    median(rhythmic_filtered$amplitude, na.rm = TRUE),
    max(rhythmic_filtered$amplitude, na.rm = TRUE),
    min(rhythmic_filtered$peak_ZT, na.rm = TRUE),
    median(rhythmic_filtered$peak_ZT, na.rm = TRUE),
    max(rhythmic_filtered$peak_ZT, na.rm = TRUE)
  )
)

write.csv(
  summary_table,
  file.path(output_dir, "rhythmic_gene_summary.csv"),
  row.names = FALSE
)


############################
# 8. Rhythmic genes by phase
############################

phase_distribution <- rhythmic_filtered %>%
  count(phase_bin, name = "number_of_genes") %>%
  complete(
    phase_bin = factor(
      phase_levels,
      levels = phase_levels
    ),
    fill = list(number_of_genes = 0)
  )

write.csv(
  phase_distribution,
  file.path(output_dir, "rhythmic_genes_by_phase.csv"),
  row.names = FALSE
)

write.csv(
  rhythmic_filtered %>%
    arrange(phase_bin, peak_ZT, padj),
  file.path(output_dir, "rhythmic_genes_final_phase_ordered.csv"),
  row.names = FALSE
)


############################
# 9. Top 10 genes per phase
############################

top_10_per_phase <- rhythmic_filtered %>%
  arrange(phase_bin, padj) %>%
  group_by(phase_bin) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top_10_per_phase,
  file.path(output_dir, "top_10_rhythmic_genes_per_phase.csv"),
  row.names = FALSE
)


############################
# 10. Phase distribution plot
############################

p_phase <- ggplot(
  rhythmic_filtered,
  aes(x = peak_ZT)
) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    closed = "left"
  ) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24)
  ) +
  labs(
    title = "Distribution of Rhythmic Gene Peak Phases",
    x = "Peak phase (ZT)",
    y = "Number of rhythmic genes"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rhythmic_gene_phase_distribution.pdf"),
  p_phase,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "rhythmic_gene_phase_distribution.png"),
  p_phase,
  width = 8,
  height = 5,
  dpi = 300
)


############################
# 11. Phase-bin bar plot
############################

p_phase_bin <- ggplot(
  phase_distribution,
  aes(
    x = phase_bin,
    y = number_of_genes
  )
) +
  geom_col() +
  labs(
    title = "Rhythmic Genes by Circadian Phase",
    x = "Peak phase",
    y = "Number of rhythmic genes"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rhythmic_genes_by_phase_bin.pdf"),
  p_phase_bin,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "rhythmic_genes_by_phase_bin.png"),
  p_phase_bin,
  width = 8,
  height = 5,
  dpi = 300
)


############################
# 12. Amplitude distribution
############################

p_amplitude <- ggplot(
  rhythmic_filtered,
  aes(x = amplitude)
) +
  geom_histogram(
    bins = 30
  ) +
  labs(
    title = "Distribution of Rhythmic Gene Amplitudes",
    x = "Amplitude",
    y = "Number of rhythmic genes"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rhythmic_gene_amplitude_distribution.pdf"),
  p_amplitude,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "rhythmic_gene_amplitude_distribution.png"),
  p_amplitude,
  width = 8,
  height = 5,
  dpi = 300
)


############################
# 13. Phase vs amplitude
############################

p_phase_amplitude <- ggplot(
  rhythmic_filtered,
  aes(
    x = peak_ZT,
    y = amplitude
  )
) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24)
  ) +
  labs(
    title = "Rhythmic Gene Phase versus Amplitude",
    x = "Peak phase (ZT)",
    y = "Amplitude"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rhythmic_phase_vs_amplitude.pdf"),
  p_phase_amplitude,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "rhythmic_phase_vs_amplitude.png"),
  p_phase_amplitude,
  width = 8,
  height = 5,
  dpi = 300
)


############################
# 14. Read VST expression
############################

cat("\n============================================\n")
cat("READING VST EXPRESSION DATA\n")
cat("============================================\n\n")

if (file.exists(vst_file)) {
  
  vst <- read.csv(
    vst_file,
    check.names = FALSE
  )
  
  # Explicitly use first column as gene ID
  colnames(vst)[1] <- "gene"
  
  vst$gene <- as.character(vst$gene)
  
  cat("VST rows:", nrow(vst), "\n")
  cat("VST columns:", ncol(vst), "\n")
  
  rhythmic_genes <- rhythmic_filtered$gene
  
  vst_rhythmic <- vst %>%
    filter(gene %in% rhythmic_genes)
  
  cat(
    "Rhythmic genes found in VST matrix:",
    nrow(vst_rhythmic),
    "\n"
  )
  
  if (nrow(vst_rhythmic) > 0) {
    
    vst_rhythmic <- vst_rhythmic %>%
      left_join(
        rhythmic_filtered %>%
          select(gene, peak_ZT, phase_bin, padj),
        by = "gene"
      ) %>%
      arrange(
        peak_ZT,
        padj
      )
    
    write.csv(
      vst_rhythmic,
      file.path(
        output_dir,
        "phase_ordered_rhythmic_gene_VST_matrix.csv"
      ),
      row.names = FALSE
    )
    
    ############################
    # VST heatmap
    ############################
    
    expression_columns <- setdiff(
      colnames(vst_rhythmic),
      c(
        "gene",
        "peak_ZT",
        "phase_bin",
        "padj"
      )
    )
    
    expression_matrix <- as.matrix(
      vst_rhythmic[, expression_columns, drop = FALSE]
    )
    
    rownames(expression_matrix) <- vst_rhythmic$gene
    
    # Convert all values to numeric
    expression_matrix <- apply(
      expression_matrix,
      2,
      as.numeric
    )
    
    rownames(expression_matrix) <- vst_rhythmic$gene
    
    ############################
    # Expected sample order
    ############################
    
    expected_ZT <- rep(
      c(0, 4, 8, 12, 16, 20),
      each = 4
    )
    
    if (ncol(expression_matrix) == 24) {
      
      sample_order <- order(expected_ZT)
      
      expression_matrix <- expression_matrix[
        ,
        sample_order,
        drop = FALSE
      ]
      
      sample_ZT <- expected_ZT[sample_order]
      
      annotation_col <- data.frame(
        ZT = factor(
          sample_ZT,
          levels = c(0, 4, 8, 12, 16, 20)
        )
      )
      
      rownames(annotation_col) <- colnames(expression_matrix)
      
      pheatmap(
        expression_matrix,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_col = annotation_col,
        show_rownames = FALSE,
        border_color = NA,
        filename = file.path(
          output_dir,
          "all_rhythmic_genes_phase_ordered_heatmap.pdf"
        ),
        width = 10,
        height = 12
      )
      
    } else {
      
      cat(
        "\nWARNING: VST matrix does not contain 24 samples.\n"
      )
      
      cat(
        "Heatmap will be generated without ZT sample annotation.\n"
      )
      
      pheatmap(
        expression_matrix,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        border_color = NA,
        filename = file.path(
          output_dir,
          "all_rhythmic_genes_phase_ordered_heatmap.pdf"
        ),
        width = 10,
        height = 12
      )
    }
    
  } else {
    
    cat(
      "No rhythmic genes were found in the VST matrix.\n"
    )
  }
  
} else {
  
  cat(
    "VST file not found. Skipping VST heatmap.\n"
  )
}


############################
# 15. Read background genes
############################

cat("\n============================================\n")
cat("READING BACKGROUND GENE SET\n")
cat("============================================\n\n")

background_genes <- NULL

if (file.exists(background_file)) {
  
  all_results <- read.csv(
    background_file,
    check.names = FALSE
  )
  
  # Explicitly use first column as gene ID
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
    "Number of genes in background:",
    length(background_genes),
    "\n"
  )
  
} else {
  
  cat(
    "Background file not found.\n"
  )
  
  cat(
    "g:Profiler will use its default background.\n"
  )
}


############################
# 16. Helper function:
#     flatten g:Profiler list columns
############################

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


############################
# 17. g:Profiler diagnostics
############################

cat("\n============================================\n")
cat("g:PROFILER SETTINGS\n")
cat("============================================\n\n")

cat(
  "Organism:",
  gprofiler_organism,
  "\n"
)

cat(
  "Sources:",
  paste(
    gprofiler_sources,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Number of query genes:",
  length(
    unique(
      rhythmic_filtered$gene
    )
  ),
  "\n"
)

if (!is.null(background_genes)) {
  
  cat(
    "Number of background genes:",
    length(background_genes),
    "\n"
  )
}


############################
# 18. Overall g:Profiler
#     enrichment
############################

cat("\n============================================\n")
cat("RUNNING OVERALL g:PROFILER ENRICHMENT\n")
cat("Species: Heterocephalus glaber\n")
cat("============================================\n\n")

rhythmic_genes <- unique(
  rhythmic_filtered$gene
)

gprofiler_success <- FALSE

gp <- tryCatch(
  
  {
    
    if (!is.null(background_genes)) {
      
      gost(
        query = rhythmic_genes,
        organism = gprofiler_organism,
        sources = gprofiler_sources,
        significant = FALSE,
        correction_method = "fdr",
        user_threshold = 0.05,
        custom_bg = background_genes,
        evcodes = TRUE
      )
      
    } else {
      
      gost(
        query = rhythmic_genes,
        organism = gprofiler_organism,
        sources = gprofiler_sources,
        significant = FALSE,
        correction_method = "fdr",
        user_threshold = 0.05,
        evcodes = TRUE
      )
    }
    
  },
  
  error = function(e) {
    
    cat(
      "\nERROR running g:Profiler:\n"
    )
    
    cat(
      conditionMessage(e),
      "\n"
    )
    
    return(NULL)
  }
)


############################
# 19. Process overall
#     enrichment
############################

if (!is.null(gp) &&
    !is.null(gp$result) &&
    nrow(gp$result) > 0) {
  
  gprofiler_success <- TRUE
  
  enrichment_table <- gp$result
  
  enrichment_table_csv <- flatten_for_csv(
    enrichment_table
  )
  
  ############################
  # Save all enrichment
  ############################
  
  write.csv(
    enrichment_table_csv,
    file.path(
      output_dir,
      "gProfiler_overall_enrichment_ALL.csv"
    ),
    row.names = FALSE
  )
  
  ############################
  # Save FDR < 0.05
  ############################
  
  significant_enrichment <- enrichment_table_csv %>%
    filter(
      p_value < 0.05
    )
  
  write.csv(
    significant_enrichment,
    file.path(
      output_dir,
      "gProfiler_overall_enrichment.csv"
    ),
    row.names = FALSE
  )
  
  
  ############################
  # GO Biological Process
  ############################
  
  go_bp <- enrichment_table_csv %>%
    filter(
      source == "GO:BP"
    )
  
  write.csv(
    go_bp,
    file.path(
      output_dir,
      "GO_Biological_Process_enrichment_all.csv"
    ),
    row.names = FALSE
  )
  
  go_bp_sig <- go_bp %>%
    filter(
      p_value < 0.05
    )
  
  write.csv(
    go_bp_sig,
    file.path(
      output_dir,
      "GO_Biological_Process_enrichment_FDR05.csv"
    ),
    row.names = FALSE
  )
  
  
  ############################
  # KEGG
  ############################
  
  kegg <- enrichment_table_csv %>%
    filter(
      source == "KEGG"
    )
  
  write.csv(
    kegg,
    file.path(
      output_dir,
      "KEGG_enrichment_all.csv"
    ),
    row.names = FALSE
  )
  
  kegg_sig <- kegg %>%
    filter(
      p_value < 0.05
    )
  
  write.csv(
    kegg_sig,
    file.path(
      output_dir,
      "KEGG_enrichment_FDR05.csv"
    ),
    row.names = FALSE
  )
  
  
  ############################
  # Diagnostics
  ############################
  
  cat("\n============================================\n")
  cat("g:PROFILER RESULTS\n")
  cat("============================================\n\n")
  
  cat(
    "Total enrichment terms returned:",
    nrow(enrichment_table_csv),
    "\n\n"
  )
  
  cat("Terms by source:\n")
  print(
    table(
      enrichment_table_csv$source
    )
  )
  
  cat("\nUnique sources returned:\n")
  print(
    unique(
      enrichment_table_csv$source
    )
  )
  
  cat(
    "\nKEGG terms returned:",
    nrow(kegg),
    "\n"
  )
  
  cat(
    "KEGG terms with FDR < 0.05:",
    nrow(kegg_sig),
    "\n"
  )
  
  cat(
    "\nGO:BP terms returned:",
    nrow(go_bp),
    "\n"
  )
  
  cat(
    "GO:BP terms with FDR < 0.05:",
    nrow(go_bp_sig),
    "\n"
  )
  
  
  ############################
  # Print top KEGG terms
  ############################
  
  if (nrow(kegg) > 0) {
    
    cat("\n============================================\n")
    cat("TOP KEGG TERMS\n")
    cat("============================================\n\n")
    
    kegg_top <- kegg %>%
      arrange(p_value) %>%
      select(
        source,
        term_id,
        term_name,
        p_value,
        intersection_size,
        intersection
      ) %>%
      head(20)
    
    print(kegg_top)
  }
  
  
  ############################
  # Overall enrichment plot
  ############################
  
  plot_data <- significant_enrichment %>%
    arrange(p_value) %>%
    head(20)
  
  if (nrow(plot_data) > 0) {
    
    p_enrichment <- ggplot(
      plot_data,
      aes(
        x = -log10(p_value),
        y = reorder(
          term_name,
          p_value
        )
      )
    ) +
      geom_point(
        size = 3
      ) +
      labs(
        title = "Overall Functional Enrichment",
        x = "-log10(FDR)",
        y = "Term"
      ) +
      theme_bw()
    
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
  
  cat(
    "\nNo g:Profiler enrichment results were returned.\n"
  )
  
  cat(
    "This may indicate a gene-ID recognition problem,\n"
  )
  
  cat(
    "an incompatible identifier type, or insufficient\n"
  )
  
  cat(
    "annotation for the submitted genes.\n"
  )
}


############################
# 20. Phase-specific
#     enrichment
############################

cat("\n============================================\n")
cat("PHASE-SPECIFIC ENRICHMENT\n")
cat("============================================\n\n")

phase_enrichment_all <- list()

for (phase in phase_levels) {
  
  phase_genes <- rhythmic_filtered %>%
    filter(
      phase_bin == phase
    ) %>%
    pull(gene) %>%
    unique()
  
  cat("\n--------------------------------------------\n")
  cat("Phase:", phase, "\n")
  cat(
    "Number of genes:",
    length(phase_genes),
    "\n"
  )
  cat("--------------------------------------------\n")
  
  if (length(phase_genes) < minimum_phase_genes) {
    
    cat(
      "Too few genes for enrichment; skipping.\n"
    )
    
    next
  }
  
  
  gp_phase <- tryCatch(
    
    {
      
      if (!is.null(background_genes)) {
        
        gost(
          query = phase_genes,
          organism = gprofiler_organism,
          sources = gprofiler_sources,
          significant = FALSE,
          correction_method = "fdr",
          user_threshold = 0.05,
          custom_bg = background_genes,
          evcodes = TRUE
        )
        
      } else {
        
        gost(
          query = phase_genes,
          organism = gprofiler_organism,
          sources = gprofiler_sources,
          significant = FALSE,
          correction_method = "fdr",
          user_threshold = 0.05,
          evcodes = TRUE
        )
      }
      
    },
    
    error = function(e) {
      
      cat(
        "ERROR for phase",
        phase,
        ":",
        conditionMessage(e),
        "\n"
      )
      
      return(NULL)
    }
  )
  
  
  if (
    is.null(gp_phase) ||
    is.null(gp_phase$result) ||
    nrow(gp_phase$result) == 0
  ) {
    
    cat(
      "No enrichment results returned for",
      phase,
      "\n"
    )
    
    next
  }
  
  
  ############################
  # Convert result
  ############################
  
  phase_table <- gp_phase$result
  
  phase_table$Phase <- phase
  
  phase_table_csv <- flatten_for_csv(
    phase_table
  )
  
  
  ############################
  # Save all phase enrichment
  ############################
  
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
  
  
  ############################
  # GO:BP
  ############################
  
  phase_go <- phase_table_csv %>%
    filter(
      source == "GO:BP"
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
  
  phase_go_sig <- phase_go %>%
    filter(
      p_value < 0.05
    )
  
  write.csv(
    phase_go_sig,
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
  
  
  ############################
  # KEGG
  ############################
  
  phase_kegg <- phase_table_csv %>%
    filter(
      source == "KEGG"
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
  
  phase_kegg_sig <- phase_kegg %>%
    filter(
      p_value < 0.05
    )
  
  write.csv(
    phase_kegg_sig,
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
  
  
  ############################
  # Diagnostics
  ############################
  
  cat(
    "GO:BP terms:",
    nrow(phase_go),
    "\n"
  )
  
  cat(
    "GO:BP FDR < 0.05:",
    nrow(phase_go_sig),
    "\n"
  )
  
  cat(
    "KEGG terms:",
    nrow(phase_kegg),
    "\n"
  )
  
  cat(
    "KEGG FDR < 0.05:",
    nrow(phase_kegg_sig),
    "\n"
  )
  
  
  ############################
  # Save for combined table
  ############################
  
  phase_enrichment_all[[phase]] <- phase_table_csv
  
  
  ############################
  # Phase enrichment plot
  ############################
  
  phase_sig <- phase_table_csv %>%
    filter(
      p_value < 0.05
    ) %>%
    arrange(
      p_value
    ) %>%
    head(20)
  
  if (nrow(phase_sig) > 0) {
    
    p_phase_enrichment <- ggplot(
      phase_sig,
      aes(
        x = -log10(p_value),
        y = reorder(
          term_name,
          p_value
        )
      )
    ) +
      geom_point(
        size = 3
      ) +
      labs(
        title = paste(
          "Functional Enrichment:",
          phase
        ),
        x = "-log10(FDR)",
        y = "Term"
      ) +
      theme_bw()
    
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
      height = 8
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
      height = 8,
      dpi = 300
    )
  }
}


############################
# 21. Combined phase
#     enrichment table
############################

if (length(phase_enrichment_all) > 0) {
  
  combined_phase_enrichment <- bind_rows(
    phase_enrichment_all
  )
  
  write.csv(
    combined_phase_enrichment,
    file.path(
      output_dir,
      "gProfiler_phase_enrichment_combined.csv"
    ),
    row.names = FALSE
  )
  
  cat("\n============================================\n")
  cat("COMBINED PHASE ENRICHMENT\n")
  cat("============================================\n\n")
  
  cat(
    "Total phase enrichment results:",
    nrow(combined_phase_enrichment),
    "\n"
  )
}


############################
# 22. Final output summary
############################

cat("\n\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n\n")

cat(
  "Species used for g:Profiler:",
  "Heterocephalus glaber",
  "\n"
)

cat(
  "g:Profiler organism ID:",
  gprofiler_organism,
  "\n"
)

cat(
  "Functional sources:",
  paste(
    gprofiler_sources,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Rhythmic genes identified:",
  nrow(rhythmic_filtered),
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n\n"
)

cat("Key outputs:\n")

cat(
  "  - rhythmic_gene_summary.csv\n"
)

cat(
  "  - rhythmic_genes_by_phase.csv\n"
)

cat(
  "  - rhythmic_genes_final_phase_ordered.csv\n"
)

cat(
  "  - top_10_rhythmic_genes_per_phase.csv\n"
)

cat(
  "  - rhythmic_gene_phase_distribution.pdf/png\n"
)

cat(
  "  - rhythmic_genes_by_phase_bin.pdf/png\n"
)

cat(
  "  - rhythmic_gene_amplitude_distribution.pdf/png\n"
)

cat(
  "  - rhythmic_phase_vs_amplitude.pdf/png\n"
)

cat(
  "  - phase_ordered_rhythmic_gene_VST_matrix.csv\n"
)

cat(
  "  - all_rhythmic_genes_phase_ordered_heatmap.pdf\n"
)

cat(
  "  - GO_Biological_Process_enrichment_all.csv\n"
)

cat(
  "  - GO_Biological_Process_enrichment_FDR05.csv\n"
)

cat(
  "  - KEGG_enrichment_all.csv\n"
)

cat(
  "  - KEGG_enrichment_FDR05.csv\n"
)

cat(
  "  - gProfiler_phase_enrichment_combined.csv\n"
)

cat("\n============================================\n")
cat("IMPORTANT INTERPRETATION NOTE\n")
cat("============================================\n\n")

cat(
  "g:Profiler p_value values are corrected using FDR.\n"
)

cat(
  "Therefore p_value < 0.05 is treated here as FDR < 0.05.\n"
)

cat(
  "The *_all.csv files contain all returned terms,\n"
)

cat(
  "including terms that are not FDR-significant.\n"
)

cat(
  "The *_FDR05.csv files contain only FDR-significant terms.\n"
)

cat(
  "\nReactome was not included because it is not an\n"
)

cat(
  "available g:Profiler source for this H. glaber analysis.\n"
)

cat(
  "\n============================================\n"
)

