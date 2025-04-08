#!/bin/bash

# Directory containing viral alignments (e.g., for herpesviruses)
ALIGN_DIR="../workdir/alignment_herpesvirus"
OUTPUT="../results/star_alignment_summary_unique_only.csv"

# Write CSV header
echo '"sample","virus","uniquely_mapped_reads","mismatch_rate_per_base"' > "$OUTPUT"

# Loop through STAR log files
for logfile in "$ALIGN_DIR"/*/Log.final.out; do
    dir=$(dirname "$logfile")
    base=$(basename "$dir")  # e.g., SRRxxxxxx__EBV
    sample=${base%%__*}
    virus=${base##*__}

    # Extract required values
    uniq_reads=$(awk '/Uniquely mapped reads number/ {print $NF}' "$logfile")
    mismatch_rate=$(awk '/Mismatch rate per base/ {print $NF}' "$logfile")

    # Handle missing values
    [[ -z "$uniq_reads" || ! "$uniq_reads" =~ ^[0-9]+$ ]] && uniq_reads=0
    [[ -z "$mismatch_rate" || "$mismatch_rate" == "nan" ]] && mismatch_rate=0

    # Output
    echo "\"$sample\",\"$virus\",\"$uniq_reads\",\"$mismatch_rate\"" >> "$OUTPUT"
done

echo "STAR alignment summary saved to: $OUTPUT"