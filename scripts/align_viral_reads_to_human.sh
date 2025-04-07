#!/bin/bash
#$ -N Align_human
#$ -cwd
#$ -o ../logs/star_align_human_$TASK_ID.out
#$ -e ../logs/star_align_human_$TASK_ID.err
#$ -pe smp 12
#$ -l h_vmem=64G
#$ -l h_rt=24:00:00
#$ -t 1-2500

source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

VIRAL_ALIGN_DIR="../workdir/alignment_herpesvirus"
HUMAN_ALIGN_DIR="../workdir/alignment_human"
VIRAL_ONLY_DIR="../workdir/viral_only"
THREADS=12
HUMAN_INDEX="../genomes/human/STAR_index"

mkdir -p "$HUMAN_ALIGN_DIR" "$VIRAL_ONLY_DIR"

dir_list=($VIRAL_ALIGN_DIR/*)

if [ $SGE_TASK_ID -gt ${#dir_list[@]} ]; then
    echo "Task ID exceeds the available directories. Exiting."
    exit 1
fi

viral_dir="${dir_list[$((SGE_TASK_ID-1))]}"
base=$(basename "$viral_dir")

bam="$viral_dir/Aligned.sortedByCoord.out.bam"
viral_log="$viral_dir/Log.final.out"

if [[ ! -f "$bam" || ! -f "$viral_log" ]]; then
    echo "Missing BAM or viral log for $base. Skipping."
    exit 1
fi

# Extract mapped viral reads counts from viral STAR log
unique_viral=$(grep "Uniquely mapped reads number" "$viral_log" | awk -F'|' '{gsub(/ /, "", $2); print $2}')
multi_viral=$(grep "Number of reads mapped to multiple loci" "$viral_log" | awk -F'|' '{gsub(/ /, "", $2); print $2}')
total_viral=$((unique_viral + multi_viral))

human_out="$HUMAN_ALIGN_DIR/$base"
mkdir -p "$human_out"

# Extract only mapped viral reads from BAM
samtools view -b -F 4 "$bam" > "$human_out/mapped_viral.bam"

# Convert mapped viral BAM to FASTQ
samtools fastq -1 "$human_out/mapped_viral_R1.fastq" -2 "$human_out/mapped_viral_R2.fastq" "$human_out/mapped_viral.bam"

# Align these reads to the human genome
STAR --runThreadN $THREADS \
     --genomeDir "$HUMAN_INDEX" \
     --readFilesIn "$human_out/mapped_viral_R1.fastq" "$human_out/mapped_viral_R2.fastq" \
     --outFileNamePrefix "$human_out/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within

human_log="$human_out/Log.final.out"
human_bam="$human_out/Aligned.sortedByCoord.out.bam"

if [[ ! -f "$human_bam" || ! -f "$human_log" ]]; then
    echo "Human alignment failed for $base."
    exit 1
fi

# Get human-aligned reads (including "too many loci")
unique_human=$(grep "Uniquely mapped reads number" "$human_log" | awk -F'|' '{gsub(/ /, "", $2); print $2}')
multi_human=$(grep "Number of reads mapped to multiple loci" "$human_log" | awk -F'|' '{gsub(/ /, "", $2); print $2}')
toomany_human=$(grep "Number of reads mapped to too many loci" "$human_log" | awk -F'|' '{gsub(/ /, "", $2); print $2}')
total_human=$((unique_human + multi_human + toomany_human))

# Viral-only calculation
viral_only_reads=$((total_viral - total_human))

viral_only_out="$VIRAL_ONLY_DIR/$base"
mkdir -p "$viral_only_out"

# Extract viral-only BAM
samtools view -b -f 4 "$human_bam" > "$viral_only_out/viral_only.bam"

# Convert viral-only BAM to FASTQ
samtools fastq -1 "$viral_only_out/viral_only_R1.fastq" -2 "$viral_only_out/viral_only_R2.fastq" "$viral_only_out/viral_only.bam"

# Write stats
cat > "$viral_only_out/stats.txt" <<EOL
Sample-Virus: $base

Viral Alignment:
  Uniquely mapped reads: $unique_viral
  Multi-mapped reads: $multi_viral
  Total viral mapped reads: $total_viral

Human Alignment:
  Uniquely mapped reads: $unique_human
  Multi-mapped reads: $multi_human
  Mapped to too many loci: $toomany_human
  Total mapped to human: $total_human

Viral-only reads: $viral_only_reads
EOL

echo "Completed $base: Viral-only results saved in $viral_only_out"