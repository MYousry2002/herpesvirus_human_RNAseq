#!/bin/bash

# Define input/output paths
ALIGN_DIR="../workdir/alignments_herpesviruses"
METADATA="../data/GSE221091/metadata.csv"
OUTPUT="../results/mapped_read_summary.csv"

# Write CSV header (all fields quoted)
echo "\"sample\",\"virus\",\"uniquely_mapped_reads\",\"multi_mapped_reads\",\"total_mapped_reads\",\"disease\",\"cell_type\"" > "$OUTPUT"

# Loop over all STAR Log.final.out files
for log in "$ALIGN_DIR"/*/Log.final.out; do
    dir=$(dirname "$log")
    base=$(basename "$dir")  # Example: SRR22753365__KSHV
    sample=${base%%__*}
    virus=${base##*__}

    # Extract STAR alignment stats
    uniq=$(grep "Uniquely mapped reads number" "$log" | cut -f2 | tr -d ' ')
    multi=$(grep "Number of reads mapped to multiple loci" "$log" | cut -f2 | tr -d ' ')
    uniq=${uniq:-0}
    multi=${multi:-0}
    total=$((uniq + multi))

    # Extract metadata using AWK with FPAT on a single line.
    IFS=$'\t' read -r meta_disease cell_type <<< "$(awk -v s="$sample" 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")" } NR > 1 && $1==s { gsub(/^"/, "", $15); gsub(/"$/, "", $15); gsub(/^"/, "", $9); gsub(/"$/, "", $9); print $15 "\t" $9; exit }' "$METADATA")"

    # Set disease based on metadata; if not "Control", use "Juvenile Dermatomyositis"
    if [ "$meta_disease" = "Control" ]; then
        disease="Control"
    else
        disease="Juvenile Dermatomyositis"
    fi

    # Write the CSV line (all fields quoted)
    echo "\"$sample\",\"$virus\",\"$uniq\",\"$multi\",\"$total\",\"$disease\",\"$cell_type\"" >> "$OUTPUT"
done

echo "Mapping summary saved to: $OUTPUT"