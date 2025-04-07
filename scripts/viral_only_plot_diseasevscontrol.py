import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file (adjust the path as needed)
df = pd.read_csv("../results/viral_only_read_summary.csv")

# Ensure the total_mapped_reads column is numeric
df['viral_only_reads'] = pd.to_numeric(df['viral_only_reads'], errors='coerce')

# Group by virus and disease, summing the total mapped reads for each group
grouped = df.groupby(['virus', 'disease'])['viral_only_reads'].sum().reset_index()

# Create a pivot table with viruses as rows and disease status as columns
pivot = grouped.pivot(index='virus', columns='disease', values='viral_only_reads').fillna(0)

# Reset index for seaborn plotting
pivot_reset = pivot.reset_index()

# Melt the pivot table so we can plot with seaborn
melted = pd.melt(pivot_reset, id_vars='virus', value_vars=pivot.columns,
                 var_name='disease', value_name='viral_only_reads')

# Set a clean style
sns.set(style="whitegrid")

# Create a grouped bar plot
plt.figure(figsize=(10, 6))
ax = sns.barplot(x="virus", y="viral_only_reads", hue="disease", data=melted)

# Add titles and labels
plt.title("Total Viral Only Mapped Reads per Virus by Disease Status")
plt.xlabel("Virus")
plt.ylabel("Total Mapped Reads")
plt.xticks(rotation=45)
plt.tight_layout()

# save the figure
plt.savefig("../results/viral_only_reads_by_virus_disease.png", dpi=300)

plt.show()