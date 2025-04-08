import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("../results/shared_unique_viruses_reads_per_sample.csv")

# Group by virus pairs and sum shared reads
summary = df.groupby(["virus_1", "virus_2"])["shared_reads"].sum().reset_index()

# Pivot to a matrix
pivot = summary.pivot(index="virus_1", columns="virus_2", values="shared_reads").fillna(0)

# Make symmetric (shared_reads is symmetric)
for v1 in pivot.index:
    for v2 in pivot.columns:
        if v1 in pivot.columns and v2 in pivot.index:
            val = (pivot.loc[v1, v2] + pivot.loc[v2, v1]) / 2
            pivot.loc[v1, v2] = val
            pivot.loc[v2, v1] = val

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(pivot, annot=True, fmt=".0f", cmap="YlOrRd", square=True, linewidths=0.5)
plt.title("Total Shared Reads Between Viruses Across All Samples")
plt.tight_layout()
plt.savefig("../results/shared_reads_heatmap.png", dpi=300)
plt.show()