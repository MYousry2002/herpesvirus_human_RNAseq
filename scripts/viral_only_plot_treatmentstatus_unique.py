import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file (adjust the path as needed)
df = pd.read_csv("../results/viral_only_mapped_read_summary_wtreatmentstatus.csv")

# Ensure the uniquely_mapped_reads column is numeric
df['uniquely_mapped_reads'] = pd.to_numeric(df['uniquely_mapped_reads'], errors='coerce')

# Group by virus and disease, summing the total mapped reads for each group
grouped = df.groupby(['virus', 'treatment_status'])['uniquely_mapped_reads'].sum().reset_index()

# Create a pivot table with viruses as rows and disease status as columns
pivot = grouped.pivot(index='virus', columns='treatment_status', values='uniquely_mapped_reads').fillna(0)

# Reset index for seaborn plotting
pivot_reset = pivot.reset_index()

# Melt the pivot table so we can plot with seaborn
melted = pd.melt(pivot_reset, id_vars='virus', value_vars=pivot.columns,
                 var_name='treatment_status', value_name='uniquely_mapped_reads')

# Set a clean style
sns.set(style="whitegrid")

# Create a grouped bar plot
plt.figure(figsize=(10, 6))
ax = sns.barplot(x="virus", y="uniquely_mapped_reads", hue="treatment_status", data=melted)

# Add titles and labels
plt.title("Uniquely Mapped Reads per Virus by Treatment Status")
plt.xlabel("Virus")
plt.ylabel("Total Mapped Reads")
plt.xticks(rotation=45)
plt.tight_layout()

# Optionally, save the figure
plt.savefig("../results/virus_only_uniquereads_treatmentstatus.png", dpi=300)

plt.show()