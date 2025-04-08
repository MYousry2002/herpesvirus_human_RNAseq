#!/bin/bash

# Define input/output paths
VIRAL_ONLY_DIR="../workdir/viral_only"
METADATA="../data/GSE221091/metadata.csv"
OUTPUT="../results/viral_only_mapped_read_summary.csv"

# Write CSV header (all fields quoted)
echo "\"sample\",\"virus\",\"uniquely_mapped_reads\",\"multi_mapped_reads\",\"total_mapped_reads\",\"disease\",\"cell_type\"" > "$OUTPUT"

# Loop over all viral_only stats files
for stats_file in "$VIRAL_ONLY_DIR"/*/stats.txt; do
    dir=$(dirname "$stats_file")
    base=$(basename "$dir")  # Example: SRR22753365__KSHV
    sample=${base%%__*}
    virus=${base##*__}

    # Extract viral-only mapping stats from stats.txt
    uniq=$(grep "Viral-only uniquely-mapped reads:" "$stats_file" | awk -F ':' '{gsub(/ /,"",$2); print $2}')
    multi=$(grep "Viral-only multi-mapped reads:" "$stats_file" | awk -F ':' '{gsub(/ /,"",$2); print $2}')
    uniq=${uniq:-0}
    multi=${multi:-0}
    total=$((uniq + multi))

    # Extract metadata (disease and cell type)
    IFS=$'\t' read -r meta_disease cell_type <<< "$(awk -v s="$sample" 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")" } NR > 1 && $1==s { gsub(/^"/, "", $15); gsub(/"$/, "", $15); gsub(/^"/, "", $9); gsub(/"$/, "", $9); print $15 "\t" $9; exit }' "$METADATA")"

    # Standardize disease labels
    if [ "$meta_disease" = "Control" ]; then
        disease="Control"
    else
        disease="Juvenile Dermatomyositis"
    fi

    # Write the CSV line
    echo "\"$sample\",\"$virus\",\"$uniq\",\"$multi\",\"$total\",\"$disease\",\"$cell_type\"" >> "$OUTPUT"
done

echo "Viral-only mapping summary saved to: $OUTPUT"