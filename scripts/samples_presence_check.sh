#!/bin/bash
# This script checks that all samples in sample_list.txt have both _1.fastq and _2.fastq files in the data directory.

SAMPLE_LIST="../data/GSE221091/sample_list.txt"
FASTQ_DIR="../data/GSE221091"
MISSING_FILE="missing_samples.txt"

# Empty the missing samples file if it exists
> "$MISSING_FILE"

while IFS= read -r sample || [ -n "$sample" ]; do
    fq1="${FASTQ_DIR}/${sample}_1.fastq"
    fq2="${FASTQ_DIR}/${sample}_2.fastq"
    
    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "$sample" >> "$MISSING_FILE"
    fi
done < "$SAMPLE_LIST"

if [ -s "$MISSING_FILE" ]; then
    echo "The following samples are missing FASTQ files (either _1 or _2):"
    cat "$MISSING_FILE"
else
    echo "All samples in $SAMPLE_LIST are present in $FASTQ_DIR."
fi