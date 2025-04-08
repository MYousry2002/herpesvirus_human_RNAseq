import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# Load dataset
df = pd.read_csv("../results/viral_only_read_summary.csv")

# Ensure numeric
df['viral_only_reads'] = pd.to_numeric(df['viral_only_reads'], errors='coerce')

# Get unique virus names
viruses = df['virus'].unique()

# Output directory
outdir = "../results/plots_by_virus"
os.makedirs(outdir, exist_ok=True)

# Loop over each virus
for virus in viruses:
    virus_df = df[df['virus'] == virus].copy()

    # Skip viruses with missing or only one disease group
    if virus_df['disease'].nunique() < 2:
        continue

    # Get groups
    group_control = virus_df[virus_df['disease'] == 'Control']['viral_only_reads']
    group_disease = virus_df[virus_df['disease'] == 'Juvenile Dermatomyositis']['viral_only_reads']

    # Skip if one of the groups is empty
    if len(group_control) == 0 or len(group_disease) == 0:
        continue

    # Mann–Whitney U test
    u_stat, p_val = mannwhitneyu(group_control, group_disease, alternative='two-sided')

    # Sample sizes
    n_control = len(group_control)
    n_disease = len(group_disease)

    # Plot
    plt.figure(figsize=(8, 6))
    order = ["Juvenile Dermatomyositis", "Control"]

    sns.boxplot(
        x="disease", y="viral_only_reads",
        data=virus_df, order=order, palette="Set2"
    )
    sns.swarmplot(
        x="disease", y="viral_only_reads",
        data=virus_df, order=order, color=".25"
    )

    # Labels
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
    plt.title(f"{virus}: Control vs. Juvenile Dermatomyositis")
    plt.xlabel("")
    plt.ylabel("Viral-Only Mapped Reads")
    plt.tight_layout()

    # Save
    plt.savefig(f"{outdir}/viral_only_{virus}_control_vs_disease.png", dpi=300)
    plt.close()