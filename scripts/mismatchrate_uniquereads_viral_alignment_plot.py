import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv("../results/mismatchrate_uniquereads_viral_alignment.csv")

# Clean up mismatch rate column
df["mismatch_rate_per_base"] = (
    df["mismatch_rate_per_base"]
    .str.replace("%", "", regex=False)
    .replace(["-nan", "nan", "NaN", ""], pd.NA)
)

# Convert to numeric (exclude bad values automatically)
df["mismatch_rate_per_base"] = pd.to_numeric(df["mismatch_rate_per_base"], errors="coerce")

# Drop rows with NA mismatch rate or 0 uniquely mapped reads
df = df.dropna(subset=["mismatch_rate_per_base"])
df = df[df["uniquely_mapped_reads"] > 0]

# Set seaborn style without grid
sns.set(style="white")

# Plot
plt.figure(figsize=(10, 6))
sns.scatterplot(
    x="uniquely_mapped_reads",
    y="mismatch_rate_per_base",
    hue="virus",
    data=df,
    palette="tab10",
    s=100,
    edgecolor="black",
    alpha=0.6  # semi-transparent dots
)

# Labels and title
plt.title("Mismatch Rate vs. Uniquely Mapped Reads")
plt.xlabel("Uniquely Mapped Reads")
plt.ylabel("Mismatch Rate per Base (%)")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

# Save and show
plt.savefig("../results/mismatchrate_uniquereads_viral_alignment_plot.png", dpi=300)
plt.show()