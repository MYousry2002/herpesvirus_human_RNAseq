import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# Load dataset
df = pd.read_csv("../results/viral_only_mapped_read_summary.csv")

# Ensure total_mapped_reads is numeric
df['total_mapped_reads'] = pd.to_numeric(df['total_mapped_reads'], errors='coerce')

# Drop rows with NA values in key columns
df = df.dropna(subset=['virus', 'disease', 'cell_type', 'total_mapped_reads'])

# Create output folder
outdir = "../results/plots_by_virus_celltype_totalreads"
os.makedirs(outdir, exist_ok=True)

# Get all unique (virus, cell_type) combinations
combinations = df[['virus', 'cell_type']].drop_duplicates()

# Loop through each combination
for _, row in combinations.iterrows():
    virus = row['virus']
    cell_type = row['cell_type']
    
    subset = df[(df['virus'] == virus) & (df['cell_type'] == cell_type)].copy()

    # Skip if <2 disease groups present
    if subset['disease'].nunique() < 2:
        continue

    # Split into groups
    group_control = subset[subset['disease'] == 'Control']['total_mapped_reads']
    group_disease = subset[subset['disease'] == 'Juvenile Dermatomyositis']['total_mapped_reads']

    # Skip if one group is empty
    if len(group_control) == 0 or len(group_disease) == 0:
        continue

    # Perform Mann–Whitney U test
    u_stat, p_val = mannwhitneyu(group_control, group_disease, alternative='two-sided')

    # Sample sizes
    n_control = len(group_control)
    n_disease = len(group_disease)

    # Plot
    plt.figure(figsize=(8, 6))
    order = ["Juvenile Dermatomyositis", "Control"]

    sns.boxplot(
        x="disease", y="total_mapped_reads",
        data=subset, order=order, palette="Set2"
    )
    sns.swarmplot(
        x="disease", y="total_mapped_reads",
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
    plt.title(f"{virus} – {cell_type}: Control vs. JDM")
    plt.xlabel("")
    plt.ylabel("Viral-Only Mapped Reads")
    plt.tight_layout()

    # Sanitize filename
    safe_celltype = cell_type.replace(" ", "_").replace("/", "_")
    safe_virus = virus.replace(" ", "_").replace("/", "_")

    # Save plot
    plt.savefig(f"{outdir}/viral_only_{safe_virus}__{safe_celltype}_control_vs_disease_totalreads.png", dpi=300)
    plt.close()