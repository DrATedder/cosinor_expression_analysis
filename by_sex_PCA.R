############################################################
# PCA PLOT GROUPED BY GONADAL SEX
############################################################

library(DESeq2)
library(ggplot2)

# 1. Set output directory
outdir <- "/home/andrew/DESeq2_rhythm"

# 2. Load existing DESeq2 object
dds <- readRDS(
  file.path(outdir, "DESeq2_24h_rhythm_dds.rds")
)

# 3. Add gonadal sex metadata
sex_info <- data.frame(
  sample = sprintf("NMR2-COL1-HYP-%02d", 1:24),
  
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

rownames(sex_info) <- sex_info$sample

# 4. Check that sample names match
if (!all(colnames(dds) %in% rownames(sex_info))) {
  stop("Some DESeq2 samples are missing from the sex metadata.")
}

# 5. Add sex to the existing DESeq2 metadata
colData(dds)$sex <- factor(
  sex_info[colnames(dds), "sex"],
  levels = c("F", "M")
)

# 6. Generate PCA data using variance-stabilised counts
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = "sex",
  returnData = TRUE
)

percentVar <- attr(pca_data, "percentVar")

# 7. Create PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = sex,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  xlab(paste0("PC1: ", round(100 * percentVar[1]), "%")) +
  ylab(paste0("PC2: ", round(100 * percentVar[2]), "%")) +
  theme_bw() +
  labs(
    title = "PCA by Gonadal Sex",
    color = "Gonadal Sex"
  )

# 8. Save plot
ggsave(
  file.path(outdir, "PCA_gonadal_sex.pdf"),
  plot = pca_plot,
  width = 8,
  height = 6
)

ggsave(
  file.path(outdir, "PCA_gonadal_sex.png"),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Sex-grouped PCA plots saved to:", outdir, "\n")
