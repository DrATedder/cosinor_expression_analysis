############################################################
# PCA: IDENTIFY GENES DRIVING PC1 AND PC2
#
# Input:
#   DESeq2_VST_expression.csv
#
# Output:
#   PCA_gene_loadings.csv
#   PCA_top_PC1_genes.csv
#   PCA_top_PC2_genes.csv
#   PCA_PC1_top_genes.pdf
#   PCA_PC2_top_genes.pdf
############################################################

library(ggplot2)

############################
# 1. Settings
############################

outdir <- "/home/andrew/DESeq2_rhythm"

input_file <- file.path(
  outdir,
  "DESeq2_VST_expression.csv"
)

TOP_N <- 20

############################
# 2. Load VST expression
############################

expr <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(expr)

storage.mode(expr) <- "numeric"

############################
# 3. Run PCA
############################

# Samples are rows; genes are columns.
# Center genes but do not scale them, matching
# the default approach used by DESeq2::plotPCA().

pca <- prcomp(
  t(expr),
  center = TRUE,
  scale. = FALSE
)

############################
# 4. Extract gene loadings
############################

loadings <- pca$rotation

gene_results <- data.frame(
  gene_id = rownames(loadings),
  PC1_loading = loadings[, 1],
  PC2_loading = loadings[, 2],
  PC1_contribution = 100 * loadings[, 1]^2,
  PC2_contribution = 100 * loadings[, 2]^2
)

############################
# 5. Save all gene loadings
############################

write.csv(
  gene_results,
  file.path(outdir, "PCA_gene_loadings.csv"),
  row.names = FALSE
)

############################
# 6. Identify top genes
############################

# Rank by absolute loading, separately for each PC.

top_PC1 <- gene_results[
  order(
    abs(gene_results$PC1_loading),
    decreasing = TRUE
  ),
]

top_PC2 <- gene_results[
  order(
    abs(gene_results$PC2_loading),
    decreasing = TRUE
  ),
]

top_PC1 <- head(top_PC1, TOP_N)

top_PC2 <- head(top_PC2, TOP_N)

############################
# 7. Save top genes
############################

write.csv(
  top_PC1,
  file.path(outdir, "PCA_top_PC1_genes.csv"),
  row.names = FALSE
)

write.csv(
  top_PC2,
  file.path(outdir, "PCA_top_PC2_genes.csv"),
  row.names = FALSE
)

############################
# 8. Plot top PC1 genes
############################

p1 <- ggplot(
  top_PC1,
  aes(
    x = reorder(gene_id, PC1_loading),
    y = PC1_loading
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top genes contributing to PC1",
    x = "Gene",
    y = "PC1 loading"
  )

ggsave(
  file.path(outdir, "PCA_PC1_top_genes.pdf"),
  p1,
  width = 8,
  height = 7
)

############################
# 9. Plot top PC2 genes
############################

p2 <- ggplot(
  top_PC2,
  aes(
    x = reorder(gene_id, PC2_loading),
    y = PC2_loading
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top genes contributing to PC2",
    x = "Gene",
    y = "PC2 loading"
  )

ggsave(
  file.path(outdir, "PCA_PC2_top_genes.pdf"),
  p2,
  width = 8,
  height = 7
)

############################
# 10. Report variance explained
############################

percent_variance <- 100 *
  pca$sdev^2 / sum(pca$sdev^2)

cat("\nPCA variance explained:\n")

cat(
  "PC1:",
  round(percent_variance[1], 2),
  "%\n"
)

cat(
  "PC2:",
  round(percent_variance[2], 2),
  "%\n"
)

cat("\nTop genes driving PC1:\n")

print(top_PC1)

cat("\nTop genes driving PC2:\n")

print(top_PC2)

cat("\nResults saved to:", outdir, "\n")

############################################################
# END
############################################################
