#!/bin/bash

ALIGN_DIR="../workdir/viral_only"
GTF_BASE="../genomes/herpesviruses"
OUTPUT="../results/viral_gene_summary.csv"

echo '"sample","virus","uniquely_mapped_reads","num_genes","top_genes"' > "$OUTPUT"

for dir in "$ALIGN_DIR"/*; do
    [[ -d "$dir" ]] || continue
    sample_virus=$(basename "$dir")
    sample=${sample_virus%%__*}
    virus=${sample_virus##*__}
    echo "Processing $sample - $virus"

    bam_file=$(find "$dir" -maxdepth 1 -name "*.bam" | head -n 1)
    [[ -f "$bam_file" ]] || { echo "❌ No BAM found for $sample_virus"; continue; }

    gtf_file="${GTF_BASE}/${virus}.gtf"
    [[ -f "$gtf_file" ]] || { echo "❌ Missing GTF: $gtf_file"; continue; }

    # Count uniquely mapped reads
    uniquely_mapped_reads=$(samtools view -c -F 0x904 "$bam_file")

    # Try featureCounts with -t gene first
    types_to_try=("gene" "CDS" "exon")
    success=0
    for ftype in "${types_to_try[@]}"; do
        featureCounts -p -a "$gtf_file" -t "$ftype" -g gene_id -o tmp.counts.txt "$bam_file" &> tmp.log
        if [[ $? -eq 0 && -f tmp.counts.txt ]]; then
            success=1
            break
        fi
    done

    if [[ $success -ne 1 ]]; then
        echo "❌ featureCounts failed on $sample_virus"
        cat tmp.log
        continue
    fi

    num_genes=$(awk 'NR > 2 && $7 > 0 {count++} END {print count+0}' tmp.counts.txt)
    top_genes=$(awk 'NR > 2 && $7 > 0 {print $1, $7}' tmp.counts.txt | sort -k2,2nr | head -n 3 | paste -sd ";" -)
    top_genes=${top_genes:-NA}

    echo "\"$sample\",\"$virus\",\"$uniquely_mapped_reads\",\"$num_genes\",\"$top_genes\"" >> "$OUTPUT"
done

rm -f tmp.counts.txt tmp.counts.txt.summary tmp.log
echo "✅ Saved viral gene summary to: $OUTPUT"