import os
import pysam
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Constants
HSV1_GENOME_LENGTH = 151974  # HSV1_KOS genome length
BIN_SIZE = 1000
NUM_BINS = (HSV1_GENOME_LENGTH + BIN_SIZE - 1) // BIN_SIZE

# Directories
VIRAL_ONLY_DIR = Path("../workdir/viral_only")

# Initialize results
binned_counts = {}

# Process HSV1_KOS samples
for bam_path in VIRAL_ONLY_DIR.glob("*__HSV1_KOS/realigned/Aligned.sortedByCoord.out.bam"):
    sample_name = bam_path.parts[-3]
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
    binned_counts[sample_name] = bin_counts

# Convert to DataFrame
bin_labels = [f"{i*BIN_SIZE}-{(i+1)*BIN_SIZE-1}" for i in range(NUM_BINS)]
df = pd.DataFrame.from_dict(binned_counts, orient="index", columns=bin_labels)

# Save
df.to_csv("../results/HSV1_KOS_read_distribution_1kb_bins.csv")

# Plot
plt.figure(figsize=(15, max(6, len(df) // 4)))
sns.heatmap(df, cmap="viridis", cbar_kws={"label": "Read Count"}, xticklabels=50)
plt.title("HSV1_KOS Read Distribution Across Genome (1kb bins)")
plt.xlabel("Genome Position (1kb bins)")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig("../results/HSV1_KOS_read_distribution_heatmap.png", dpi=300)
plt.show()