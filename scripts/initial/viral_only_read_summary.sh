#!/bin/bash

# Define input/output paths
VIRAL_ONLY_DIR="../workdir/viral_only"
METADATA="../data/GSE221091/metadata.csv"
OUTPUT="../results/viral_only_read_summary.csv"

# Write CSV header (all fields quoted)
echo '"sample","virus","viral_only_reads","disease","cell_type"' > "$OUTPUT"

# Loop over all viral-only stats.txt files
for stats in "$VIRAL_ONLY_DIR"/*/stats.txt; do
    dir=$(dirname "$stats")
    base=$(basename "$dir")  # Example: SRR22753365__KSHV
    sample=${base%%__*}
    virus=${base##*__}

    # Extract viral-only reads from stats.txt
    viral_only_reads=$(grep "Viral-only reads:" "$stats" | awk '{print $3}')
    viral_only_reads=${viral_only_reads:-0}

    # Extract metadata using AWK with FPAT
    IFS=$'\t' read -r meta_disease cell_type <<< "$(awk -v s="$sample" 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")" } NR > 1 && $1==s { gsub(/^"/, "", $15); gsub(/"$/, "", $15); gsub(/^"/, "", $9); gsub(/"$/, "", $9); print $15 "\t" $9; exit }' "$METADATA")"

    # Set disease based on metadata; if not "Control", use "Juvenile Dermatomyositis"
    if [ "$meta_disease" = "Control" ]; then
        disease="Control"
    else
        disease="Juvenile Dermatomyositis"
    fi

    # Write the CSV line (all fields quoted)
    echo "\"$sample\",\"$virus\",\"$viral_only_reads\",\"$disease\",\"$cell_type\"" >> "$OUTPUT"
done

echo "Viral-only mapping summary saved to: $OUTPUT"
