#!/bin/bash

# Define paths
VIRAL_ONLY_DIR="../workdir/viral_only"
METADATA="../data/GSE221091/metadata.csv"
PATIENT_METADATA="../data/GSE221091/patient_metadata.csv"
OUTPUT="../results/viral_only_read_summary_treatment_comparison.csv"

# Write header
echo '"sample","virus","viral_only_reads","disease","treatment_status"' > "$OUTPUT"

# Loop through each viral-only stats file
for stats in "$VIRAL_ONLY_DIR"/*/stats.txt; do
    dir=$(dirname "$stats")
    base=$(basename "$dir")           # e.g., SRR22753365__HSV1_KOS
    sample=${base%%__*}               # SRR ID
    virus=${base##*__}                # Virus name

    # Get viral-only read count
    viral_only_reads=$(grep "Viral-only reads:" "$stats" | awk '{print $3}')
    viral_only_reads=${viral_only_reads:-0}

    # Get GSM ID from metadata.csv (column 27)
    gsm_id=$(awk -v s="$sample" 'BEGIN{FPAT="([^,]+)|(\"[^\"]+\")"} NR>1 && $1==s {print $27; exit}' "$METADATA")

    # Extract relevant metadata from patient_metadata.csv using GSM ID
    IFS=',' read -r id name condition treatment_status cell_type <<< "$(awk -F',' -v gsm="$gsm_id" 'NR > 1 && $1 == gsm {print $0; exit}' "$PATIENT_METADATA")"

    # Only include JDM samples with valid treatment status
    if [[ "$condition" != "Juvenile Dermatomyositis" || "$treatment_status" == "null" || -z "$treatment_status" ]]; then
        continue
    fi

    # Write to output CSV
    echo "\"$sample\",\"$virus\",\"$viral_only_reads\",\"$condition\",\"$treatment_status\"" >> "$OUTPUT"
done

echo "✔️ Treatment-based viral-only summary saved to: $OUTPUT"