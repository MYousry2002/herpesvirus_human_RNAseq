import os
import pysam
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Constants
HSV1_GENOME_LENGTH = 151974
BIN_SIZE = 1000
NUM_BINS = (HSV1_GENOME_LENGTH + BIN_SIZE - 1) // BIN_SIZE

# Paths
VIRAL_ONLY_DIR = Path("../workdir/viral_only")
METADATA_PATH = "../data/GSE221091/metadata.csv"
OUT_CSV = "../results/HSV1_KOS_read_distribution_1kb_bins.csv"
OUT_PNG = "../results/HSV1_KOS_read_distribution_heatmap_genename_groupedby_disease_celltype.png"

# Load metadata
meta = pd.read_csv(METADATA_PATH)
meta.columns = meta.columns.str.strip()
sample_to_disease = dict(zip(meta["Run"], meta["disease"]))
sample_to_celltype = dict(zip(meta["Run"], meta["cell_type"]))

# Gather data
binned_counts = []
multi_index = []

for bam_path in VIRAL_ONLY_DIR.glob("*__HSV1_KOS/realigned/Aligned.sortedByCoord.out.bam"):
    full_sample_name = bam_path.parts[-3]
    sample_id = full_sample_name.split("__")[0]

    disease = sample_to_disease.get(sample_id, "Unknown")
    celltype = sample_to_celltype.get(sample_id, "Unknown")

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    bin_counts = np.zeros(NUM_BINS, dtype=int)

    for read in bamfile.fetch():
        if not read.is_unmapped:
            bin_index = read.reference_start // BIN_SIZE
            if bin_index < NUM_BINS:
                bin_counts[bin_index] += 1

    bamfile.close()
    binned_counts.append(bin_counts)
    multi_index.append((disease, celltype, full_sample_name))

# Build DataFrame
df = pd.DataFrame(
    binned_counts,
    index=pd.MultiIndex.from_tuples(multi_index, names=["Disease", "CellType", "Sample"]),
    columns=[f"{i*BIN_SIZE}-{(i+1)*BIN_SIZE-1}" for i in range(NUM_BINS)]
)
df.sort_index(inplace=True)
df.to_csv(OUT_CSV)

# Plotting
fig, ax = plt.subplots(figsize=(15, max(6, len(df) // 4)))
sns.set(style="white")

# Create the heatmap
heatmap = sns.heatmap(
    df,
    cmap="viridis",
    cbar_kws={"label": "Read Count"},
    xticklabels=1,
    ax=ax
)

# Move colorbar to the far right
cbar = heatmap.collections[0].colorbar
cbar_ax = cbar.ax
fig.subplots_adjust(left=0.1, right=0.86, top=0.93, bottom=0.06)
cbar_ax.set_position([0.88, 0.1, 0.015, 0.8])  # [left, bottom, width, height]

# X-axis ticks every 1kb, label every 10kb
tick_positions = np.arange(0, NUM_BINS, 1)
tick_labels = [f"{i}kb" if i % 10 == 0 else "" for i in tick_positions]
ax.set_xticks(tick_positions + 0.5)
ax.set_xticklabels(tick_labels, rotation=90, fontsize=8)

# Y-axis: sample names
ax.set_yticks(np.arange(len(df)) + 0.5)
ax.set_yticklabels(df.index.get_level_values("Sample"), fontsize=6)

# Draw horizontal lines to separate disease and cell types
disease_labels = df.index.get_level_values("Disease")
celltype_labels = df.index.get_level_values("CellType")

# Find boundaries where disease changes
disease_boundaries = np.where(disease_labels[:-1] != disease_labels[1:])[0] + 1

# Find boundaries where cell type changes (within the same disease group)
celltype_boundaries = [
    i + 1 for i in range(len(celltype_labels) - 1)
    if celltype_labels[i] != celltype_labels[i + 1] and disease_labels[i] == disease_labels[i + 1]
]

# Draw main (prominent) horizontal lines for disease separation
for b in disease_boundaries:
    ax.axhline(b, color="white", linestyle="-", linewidth=2.2)

# Draw lighter horizontal lines for cell type separation
for b in celltype_boundaries:
    ax.axhline(b, color="white", linestyle="--", linewidth=1)

# Add group annotations (disease | cell_type)
for i, (disease, celltype, _) in enumerate(df.index):
    ax.text(
        len(df.columns) + 1.5, i + 0.5,
        f"{disease} | {celltype}",
        va="center", fontsize=6, color="black"
    )

# final touches
ax.set_title("HSV1_KOS Read Distribution (1kb bins)\nGrouped by Disease and Cell Type", pad=50)
ax.set_xlabel("Genome Position", labelpad=8)
ax.set_ylabel("Sample", labelpad=8)

# === Add annotated gene regions to the heatmap ===
GTF_PATH = "../genomes/herpesviruses/HSV1_KOS.gtf"  # adjust path as needed

gene_regions = []
with open(GTF_PATH) as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attributes = parts
        if feature != "gene":
            continue
        start = int(start)
        end = int(end)
        bin_start = start // BIN_SIZE
        bin_end = end // BIN_SIZE

        # Try to extract gene name
        gene_name = "unknown"
        for field in attributes.strip().split(";"):
            if "gene_name" in field or "Name" in field:
                gene_name = field.split('"')[1]
                break

        gene_regions.append((bin_start, bin_end, gene_name))

# === Plot gene regions and label ALL gene names ===
for bin_start, bin_end, gene_name in gene_regions:
    # Shaded gene region
    ax.axvspan(bin_start, bin_end, color="white", alpha=0.12)

    # Label gene name
    ax.text(
        (bin_start + bin_end) / 2, -1,  # Just above the heatmap
        gene_name,
        ha="center", va="bottom",
        fontsize=5, rotation=90, color="black", clip_on=False
    )

# Save
plt.savefig(OUT_PNG, dpi=300)
plt.show()