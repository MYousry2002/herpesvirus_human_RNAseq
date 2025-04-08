#!/bin/bash
#$ -N Rerun_combined_pipeline
#$ -cwd
#$ -o ../logs/rerun_combined_pipeline_$TASK_ID.out
#$ -e ../logs/rerun_combined_pipeline_$TASK_ID.err
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
#$ -t 1-250

# Load conda environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Directories
VIRAL_ALIGN_DIR="../workdir/alignment_herpesvirus"
HUMAN_ALIGN_DIR="../workdir/alignment_human"
VIRAL_ONLY_DIR="../workdir/viral_only"
VIRAL_INDEX_BASE="../genomes/herpesviruses/STAR_indices"
HUMAN_INDEX="../genomes/human/STAR_index"
THREADS=4

# Get only combined genome alignments
combined_dirs=($VIRAL_ALIGN_DIR/*__combined_herpes_viruses)
TOTAL_TASKS=${#combined_dirs[@]}

# Validate task ID
if [ -z "$SGE_TASK_ID" ] || [ "$SGE_TASK_ID" -gt "$TOTAL_TASKS" ]; then
    echo "Invalid SGE_TASK_ID: $SGE_TASK_ID"
    exit 1
fi

# Define this task's sample
viral_dir="${combined_dirs[$((SGE_TASK_ID-1))]}"
base=$(basename "$viral_dir")
sample=${base%%__*}
virus=${base##*__}

bam="$viral_dir/Aligned.sortedByCoord.out.bam"
viral_log="$viral_dir/Log.final.out"
human_out="$HUMAN_ALIGN_DIR/$base"
viral_only_out="$VIRAL_ONLY_DIR/$base"
realigned_out="$viral_only_out/realigned"

# Clean any previous output
rm -rf "$human_out" "$viral_only_out"
mkdir -p "$human_out" "$realigned_out"

# Check files exist
if [[ ! -f "$bam" || ! -f "$viral_log" ]]; then
    echo "Missing BAM or STAR log for $base"
    exit 1
fi

# --- Viral STAR stats ---
unique_viral=$(grep "Uniquely mapped reads number" "$viral_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
multi_viral=$(grep "Number of reads mapped to multiple loci" "$viral_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
total_viral=$((unique_viral + multi_viral))

# --- Extract mapped reads from viral BAM ---
samtools view -b -F 4 "$bam" > "$human_out/mapped_viral.bam"
samtools sort -n "$human_out/mapped_viral.bam" -o "$human_out/mapped_viral.query.bam"
samtools fastq -1 "$human_out/mapped_viral_R1.fastq" -2 "$human_out/mapped_viral_R2.fastq" -0 /dev/null -s /dev/null -n "$human_out/mapped_viral.query.bam"

# --- Align to human genome ---
STAR --runThreadN "$THREADS" \
     --genomeDir "$HUMAN_INDEX" \
     --readFilesIn "$human_out/mapped_viral_R1.fastq" "$human_out/mapped_viral_R2.fastq" \
     --outFileNamePrefix "$human_out/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 9999

human_log="$human_out/Log.final.out"
human_bam="$human_out/Aligned.sortedByCoord.out.bam"

# Validate human alignment
if [[ ! -f "$human_log" || ! -f "$human_bam" ]]; then
    echo "Human alignment failed for $base"
    exit 1
fi

# --- Human STAR stats ---
unique_human=$(grep "Uniquely mapped reads number" "$human_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
multi_human=$(grep "Number of reads mapped to multiple loci" "$human_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
total_human=$((unique_human + multi_human))

# --- Extract human-mapped read names ---
samtools view -F 4 "$human_bam" | cut -f1 | sort -u | sed 's/^/@/' > "$viral_only_out/human_mapped_readnames.txt"

# --- Filter unmapped reads from FASTQs ---
filter_fastq() {
  input_fastq="$1"
  output_fastq="$2"
  grep_list="$3"

  awk -v list="$grep_list" '
    BEGIN {
      while ((getline line < list) > 0) ids[line] = 1
    }
    {
      header = $0
      getline seq
      getline plus
      getline qual
      if (!(header in ids)) {
        print header; print seq; print plus; print qual
      }
    }
  ' "$input_fastq" > "$output_fastq"
}

filter_fastq "$human_out/mapped_viral_R1.fastq" "$viral_only_out/viral_only_R1.fastq" "$viral_only_out/human_mapped_readnames.txt"
filter_fastq "$human_out/mapped_viral_R2.fastq" "$viral_only_out/viral_only_R2.fastq" "$viral_only_out/human_mapped_readnames.txt"

# --- Re-align viral-only reads to combined genome ---
STAR --runThreadN "$THREADS" \
     --genomeDir "$VIRAL_INDEX_BASE/combined_herpes_viruses" \
     --readFilesIn "$viral_only_out/viral_only_R1.fastq" "$viral_only_out/viral_only_R2.fastq" \
     --outFileNamePrefix "$realigned_out/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMismatchNoverReadLmax 0.1 \
     --outFilterMismatchNmax 20 \
     --outFilterMultimapNmax 9999
     

# --- Re-alignment stats ---
realigned_log="$realigned_out/Log.final.out"
if [[ -f "$realigned_log" ]]; then
    viral_only_unique=$(grep "Uniquely mapped reads number" "$realigned_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
    viral_only_multi=$(grep "Number of reads mapped to multiple loci" "$realigned_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
    viral_only_total=$((viral_only_unique + viral_only_multi))
else
    viral_only_unique=0
    viral_only_multi=0
    viral_only_total=0
fi

# --- Save summary stats ---
cat > "$viral_only_out/stats.txt" <<EOL
Sample-Virus: $base

Viral Alignment:
  Uniquely mapped reads:  $unique_viral
  Multi-mapped reads:     $multi_viral
  Total viral mapped reads: $total_viral

Human Alignment:
  Uniquely mapped reads:  $unique_human
  Multi-mapped reads:     $multi_human
  Total mapped to human:  $total_human

Viral-only uniquely-mapped reads: $viral_only_unique
Viral-only multi-mapped reads:    $viral_only_multi
Viral-only total reads:           $viral_only_total
EOL

echo "✅ Completed: $base"