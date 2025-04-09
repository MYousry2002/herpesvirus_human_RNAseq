import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# Load dataset
df = pd.read_csv("../results/viral_only_mapped_read_summary.csv")

# Ensure uniquely mapped reads column is numeric
df['uniquely_mapped_reads'] = pd.to_numeric(df['uniquely_mapped_reads'], errors='coerce')

# Drop rows with NA values in important columns
df = df.dropna(subset=['virus', 'disease', 'cell_type', 'uniquely_mapped_reads'])

# Output directory
outdir = "../results/plots_by_virus_celltype_uniquereads"
os.makedirs(outdir, exist_ok=True)

# Get unique virus–cell type combinations
combinations = df[['virus', 'cell_type']].drop_duplicates()

# Loop through each virus–cell type combo
for _, row in combinations.iterrows():
    virus = row['virus']
    cell_type = row['cell_type']
    
    subset = df[(df['virus'] == virus) & (df['cell_type'] == cell_type)].copy()

    if subset['disease'].nunique() < 2:
        continue

    group_control = subset[subset['disease'] == 'Control']['uniquely_mapped_reads']
    group_disease = subset[subset['disease'] == 'Juvenile Dermatomyositis']['uniquely_mapped_reads']

    if len(group_control) == 0 or len(group_disease) == 0:
        continue

    u_stat, p_val = mannwhitneyu(group_control, group_disease, alternative='two-sided')

    n_control = len(group_control)
    n_disease = len(group_disease)

    # Plot
    plt.figure(figsize=(8, 6))
    order = ["Juvenile Dermatomyositis", "Control"]

    sns.boxplot(
        x="disease", y="uniquely_mapped_reads",
        data=subset, order=order, palette="Set2"
    )
    sns.swarmplot(
        x="disease", y="uniquely_mapped_reads",
        data=subset, order=order, color=".25"
    )

    plt.xticks(
        [0, 1],
        [
            f"Juvenile Dermatomyositis (n={n_disease})",
            f"Control (n={n_control})"
        ]
    )
    plt.text(
        0.5, 0.9, f"Mann–Whitney p={p_val:.3e}",
        transform=plt.gca().transAxes, ha='center', va='center'
    )
    plt.title(f"{virus} – {cell_type}: Unique Reads (Control vs. JDM)")
    plt.xlabel("")
    plt.ylabel("Viral-Only Uniquely Mapped Reads")
    plt.tight_layout()

    safe_celltype = cell_type.replace(" ", "_").replace("/", "_")
    safe_virus = virus.replace(" ", "_").replace("/", "_")

    plt.savefig(f"{outdir}/viral_only_{safe_virus}__{safe_celltype}_control_vs_disease_uniquereads.png", dpi=300)
    plt.close()