import os
from collections import defaultdict
from pathlib import Path
import csv

VIRAL_ONLY_DIR = "../workdir/viral_only"
OUTPUT_CSV = "../results/shared_unique_viruses_reads_per_sample.csv"

# Group paths by sample
sample_virus_fastqs = defaultdict(dict)

for path in Path(VIRAL_ONLY_DIR).glob("*__*/viral_only_R1.fastq"):
    base = path.parent.name
    sample, virus = base.split("__", 1)
    sample_virus_fastqs[sample][virus] = path

# Extract read names and compare
results = []
all_fields = set(["sample", "virus_1", "virus_2", "shared_reads"])

for sample, virus_fastqs in sample_virus_fastqs.items():
    virus_reads = {}
    for virus, fq_path in virus_fastqs.items():
        with open(fq_path, "r") as f:
            reads = set(line.strip() for i, line in enumerate(f) if i % 4 == 0)
            virus_reads[virus] = reads

    viruses = list(virus_reads.keys())
    for i, v1 in enumerate(viruses):
        for j in range(i + 1, len(viruses)):
            v2 = viruses[j]
            reads1 = virus_reads[v1]
            reads2 = virus_reads[v2]

            shared = reads1 & reads2
            only_v1 = reads1 - reads2
            only_v2 = reads2 - reads1

            results.append({
                "sample": sample,
                "virus_1": v1,
                "virus_2": v2,
                "shared_reads": len(shared),
                f"unique_to_{v1}": len(only_v1),
                f"unique_to_{v2}": len(only_v2),
                f"total_reads_{v1}": len(reads1),
                f"total_reads_{v2}": len(reads2),
            })

    # Add self-comparisons
    for virus, reads in virus_reads.items():
        results.append({
            "sample": sample,
            "virus_1": virus,
            "virus_2": virus,
            "shared_reads": len(reads),
            f"unique_to_{virus}": 0,
            f"total_reads_{virus}": len(reads)
        })

# Collect all field names used across all rows
for row in results:
    all_fields.update(row.keys())

all_fields = sorted(all_fields)

# Write to CSV
with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_fields)
    writer.writeheader()
    for row in results:
        writer.writerow(row)
        
print(f"Shared/unique read comparison saved to: {OUTPUT_CSV}")