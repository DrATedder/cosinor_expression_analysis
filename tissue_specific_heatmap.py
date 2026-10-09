import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load expression data
df = pd.read_csv(input_file)

# Remove rows without gene names
df = df.dropna(subset=["Genes"])
df["Genes"] = df["Genes"].astype(str).str.strip()
df = df[df["Genes"] != ""]

# Set gene names as the row index
df = df.set_index("Genes")

# Convert expression values to numeric
df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

# Log-transform counts
log_df = np.log2(df + 1)

# Define gene groups
hypothalamic_genes = [
    "Agrp", "Npy", "Pomc", "Sst", "Sim1", "Avp", "Oxt",
    "Pdyn", "Ghrh", "Crh", "Pmch", "Kiss1", "Trh",
    "Otp", "Hcrt", "Cartpt"
]

mouse_model_genes = [
    "Fezf1", "Gal", "Gabrq", "Slc18a2", "Magel2",
    "Slc6a3", "Gpx3", "Ngb", "Baiap3"
]

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

# Order genes by group, retaining only genes present in the data
ordered_genes = []

for group in [
    hypothalamic_genes,
    mouse_model_genes,
    negative_marker_genes,
    nearby_brain_genes
]:
    ordered_genes.extend(
        gene for gene in group
        if gene in log_df.index and gene not in ordered_genes
    )

# Add any remaining genes
ordered_genes.extend(
    gene for gene in log_df.index
    if gene not in ordered_genes
)

log_df = log_df.loc[ordered_genes]

# Create heatmap with original colours and white background
fig, ax = plt.subplots(figsize=(20, 12))

im = ax.imshow(
    log_df.values,
    aspect="auto",
    interpolation="nearest",
    cmap="viridis"
)

# Axis labels
ax.set_yticks(np.arange(len(log_df.index)))
ax.set_yticklabels(log_df.index, fontsize=11)

ax.set_xticks(np.arange(len(log_df.columns)))
ax.set_xticklabels(
    [f"{i:02d}" for i in range(1, len(log_df.columns) + 1)],
    fontsize=10
)

ax.set_xlabel("Sample", fontsize=12)
ax.set_ylabel("Gene", fontsize=12)
ax.set_title(
    "Hypothalamic and Mouse-Model Gene Expression",
    fontsize=15,
    fontweight="bold",
    pad=15
)

# Separate gene groups
n_hyp = sum(g in hypothalamic_genes for g in log_df.index)
n_mouse = sum(g in mouse_model_genes for g in log_df.index)

if 0 < n_hyp < len(log_df):
    ax.axhline(n_hyp - 0.5, color="white", linewidth=1.5)

if 0 < n_hyp + n_mouse < len(log_df):
    ax.axhline(n_hyp + n_mouse - 0.5, color="white", linewidth=1.5)

# Colour scale
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.025)
cbar.set_label("log₂(count + 1)", fontsize=11)

# Standard white figure background and black labels
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Layout
fig.tight_layout()

# Save and display
plt.savefig(
    "hypothalamus_mouse_model_heatmap.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)

plt.savefig(
    "hypothalamus_mouse_model_heatmap.pdf",
    bbox_inches="tight",
    facecolor="white"
)

plt.show()
