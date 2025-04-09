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
OUT_PNG = "../results/HSV1_KOS_read_distribution_heatmap_groupedbydisease.png"

# Load metadata
meta = pd.read_csv(METADATA_PATH)
meta.columns = meta.columns.str.strip()
sample_to_disease = dict(zip(meta["Run"], meta["disease"]))

# Read binning
grouped_counts = {"Juvenile Dermatomyositis": {}, "Control": {}}

for bam_path in VIRAL_ONLY_DIR.glob("*__HSV1_KOS/realigned/Aligned.sortedByCoord.out.bam"):
    full_sample_name = bam_path.parts[-3]
    sample_id = full_sample_name.split("__")[0]
    disease = sample_to_disease.get(sample_id, "Unknown")

    if disease not in grouped_counts:
        continue

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    bin_counts = np.zeros(NUM_BINS, dtype=int)

    for read in bamfile.fetch():
        if read.is_unmapped:
            continue
        pos = read.reference_start
        bin_index = pos // BIN_SIZE
        if bin_index < NUM_BINS:
            bin_counts[bin_index] += 1

    bamfile.close()
    grouped_counts[disease][f"{full_sample_name}"] = bin_counts

# Build dataframes with multi-index for group labels
df_jdm = pd.DataFrame.from_dict(grouped_counts["Juvenile Dermatomyositis"], orient="index")
df_control = pd.DataFrame.from_dict(grouped_counts["Control"], orient="index")

# Add multi-index for clearer group separation
df_jdm.index = pd.MultiIndex.from_product([["JDM"], df_jdm.index])
df_control.index = pd.MultiIndex.from_product([["Control"], df_control.index])

# Combine and rename columns
df = pd.concat([df_jdm, df_control])
df.columns = [f"{i*BIN_SIZE}-{(i+1)*BIN_SIZE-1}" for i in range(NUM_BINS)]
df.to_csv(OUT_CSV)

# Plot
plt.figure(figsize=(15, max(6, len(df) // 4)))
ax = sns.heatmap(df, cmap="viridis", cbar_kws={"label": "Read Count"}, xticklabels=1)

# Set xticks every 10kb
tick_positions = np.arange(0, NUM_BINS, 1)
tick_labels = [f"{i}kb" if i % 10 == 0 else "" for i in tick_positions]
plt.xticks(tick_positions + 0.5, tick_labels, rotation=90, fontsize=8)

# Replace y-tick labels with sample names only
ytick_labels = [idx[1] for idx in df.index]
plt.yticks(np.arange(len(ytick_labels)) + 0.5, ytick_labels, fontsize=6)

# Add horizontal line between groups
plt.axhline(len(df_jdm), color="white", linestyle="--", linewidth=1.2)

# Add group labels to the right of the plot
for i, (group, _) in enumerate(df.index):
    y = i + 0.5
    if group == "JDM":
        plt.text(len(df.columns) + 1, y, "JDM", va="center", fontsize=6, color="black")
    else:
        plt.text(len(df.columns) + 1, y, "Control", va="center", fontsize=6, color="black")

# Final touches
plt.title("HSV1_KOS Read Distribution Across Genome (1kb bins) by Disease")
plt.xlabel("Genome Position")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
plt.show()