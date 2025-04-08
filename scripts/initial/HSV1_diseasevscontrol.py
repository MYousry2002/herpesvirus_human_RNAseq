import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# 1. Load and filter data
df = pd.read_csv("../results/mapped_read_summary.csv")
hsv1_df = df[df['virus'] == 'HSV1_KOS'].copy()

# Ensure total_mapped_reads is numeric
hsv1_df['total_mapped_reads'] = pd.to_numeric(hsv1_df['total_mapped_reads'], errors='coerce')

# 2. Separate groups
group_control = hsv1_df[hsv1_df['disease'] == 'Control']['total_mapped_reads']
group_disease = hsv1_df[hsv1_df['disease'] == 'Juvenile Dermatomyositis']['total_mapped_reads']

# 3. Mann–Whitney U test
u_stat, p_val = mannwhitneyu(group_control, group_disease, alternative='two-sided')

# Count sample sizes
n_control = len(group_control)
n_disease = len(group_disease)

print("Control sample size:", n_control)
print("Juvenile Dermatomyositis sample size:", n_disease)
print("Mann–Whitney U statistic:", u_stat)
print("p-value:", p_val)

# 4. Plot (boxplot + swarmplot)
plt.figure(figsize=(8, 6))

# Force the order: JD first, then Control
order = ["Juvenile Dermatomyositis", "Control"]

sns.boxplot(
    x="disease", y="total_mapped_reads",
    data=hsv1_df, order=order, palette="Set2"
)
sns.swarmplot(
    x="disease", y="total_mapped_reads",
    data=hsv1_df, order=order, color=".25"
)

# Replace x-tick labels with disease + sample size
plt.xticks(
    [0, 1],
    [
        f"Juvenile Dermatomyositis (n={n_disease})",
        f"Control (n={n_control})"
    ]
)

# Place the p-value text inside the plot
plt.text(
    0.5, 0.9, f"Mann–Whitney p={p_val:.3e}",
    transform=plt.gca().transAxes, ha='center', va='center'
)

# Final formatting
plt.title("HSV1_KOS: Control vs. Juvenile Dermatomyositis")
plt.xlabel("")
plt.ylabel("Total Mapped Reads")
plt.tight_layout()

# Save and show
plt.savefig("../results/HSV1_KOS_control_vs_disease.png", dpi=300)
plt.show()