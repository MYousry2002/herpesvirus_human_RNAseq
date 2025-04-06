#!/bin/bash
#$ -N download_data
#$ -cwd
#$ -o ../logs/download_data$TASK_ID.out
#$ -e ../logs/download_data$TASK_ID.err
#$ -pe smp 12
#$ -l h_vmem=64G
#$ -l h_rt=08:00:00
#$ -t 1-50

# Load conda env and tools
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Define paths
SRR_LIST="../data/GSE221091/sample_list.txt"
OUTDIR="../data/GSE221091"
mkdir -p "$OUTDIR"

Prefetch all SRRs
echo "Prefetching all SRA files ..."
cat "$SRR_LIST" | while read SRR; do
    echo "Prefetching $SRR ..."
    prefetch "$SRR"
done

# Convert to FASTQ using fasterq-dump
echo "Converting SRA to FASTQ ..."
cat "$SRR_LIST" | while read SRR; do
    echo "Running fasterq-dump for $SRR ..."
    fasterq-dump "$SRR" -O "$OUTDIR" --split-files --threads 6
done

echo "All FASTQ files downloaded to $OUTDIR"