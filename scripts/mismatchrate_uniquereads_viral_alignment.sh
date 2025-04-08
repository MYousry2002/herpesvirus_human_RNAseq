#!/bin/bash

# Directory containing viral alignments (e.g., for herpesviruses)
ALIGN_DIR="../workdir/viral_only"
OUTPUT="../results/mismatchrate_uniquereads_viral_alignment.csv"

# Write CSV header
echo '"sample","virus","uniquely_mapped_reads","mismatch_rate_per_base"' > "$OUTPUT"

# Loop through STAR log files
for logfile in "$ALIGN_DIR"/*/realigned/Log.final.out; do
    # Skip if no matching files
    [[ ! -f "$logfile" ]] && continue

    dir=$(dirname "$logfile")
    parent_dir=$(basename "$(dirname "$dir")")  # gets SRRxxxx__Virus
    sample=${parent_dir%%__*}
    virus=${parent_dir##*__}

    # Extract values
    uniq_reads=$(awk '/Uniquely mapped reads number/ {print $NF}' "$logfile")
    mismatch_rate=$(awk '/Mismatch rate per base/ {print $NF}' "$logfile")

    # Handle missing values
    [[ -z "$uniq_reads" || ! "$uniq_reads" =~ ^[0-9]+$ ]] && uniq_reads=0
    [[ -z "$mismatch_rate" || "$mismatch_rate" == "nan" ]] && mismatch_rate=0

    echo "\"$sample\",\"$virus\",\"$uniq_reads\",\"$mismatch_rate\"" >> "$OUTPUT"
done

echo "STAR alignment summary saved to: $OUTPUT"