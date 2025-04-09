import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
df = pd.read_csv("../results/viral_only_mapped_read_summary.csv")

# Ensure the total_mapped_reads column is numeric
df['total_mapped_reads'] = pd.to_numeric(df['total_mapped_reads'], errors='coerce')

# Group by virus and cell_type (lowercase)
grouped = df.groupby(['virus', 'cell_type'])['total_mapped_reads'].sum().reset_index()

# Pivot for plotting
pivot = grouped.pivot(index='virus', columns='cell_type', values='total_mapped_reads').fillna(0)

# Reset index for seaborn
pivot_reset = pivot.reset_index()

# Melt for seaborn barplot
melted = pd.melt(pivot_reset, id_vars='virus', value_vars=pivot.columns,
                 var_name='cell_type', value_name='total_mapped_reads')

# Set style
sns.set(style="whitegrid")

# Plot
plt.figure(figsize=(12, 6))
ax = sns.barplot(x="virus", y="total_mapped_reads", hue="cell_type", data=melted)

# Labels
plt.title("Total Mapped Reads per Virus by Cell Type")
plt.xlabel("Virus")
plt.ylabel("Total Mapped Reads")
plt.xticks(rotation=45)
plt.tight_layout()

# Save
plt.savefig("../results/virus_only_totalreads_celltype.png", dpi=300)

plt.show()