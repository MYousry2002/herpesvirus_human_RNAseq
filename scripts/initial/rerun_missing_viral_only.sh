#!/bin/bash
#$ -N Align_human_virus
#$ -cwd
#$ -o ../logs/ret_star_align_human_herpesvirus_$TASK_ID.out
#$ -e ../logs/ret_star_align_human_herpesvirus_$TASK_ID.err
#$ -pe smp 12
#$ -l h_vmem=64G
#$ -l h_rt=12:00:00
#$ -t 1-200

source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Paths
VIRAL_ALIGN_DIR="../workdir/alignment_herpesvirus"
HUMAN_ALIGN_DIR="../workdir/alignment_human"
VIRAL_ONLY_DIR="../workdir/viral_only"
VIRAL_INDEX_BASE="../genomes/herpesviruses/STAR_indices"
HUMAN_INDEX="../genomes/human/STAR_index"
MISSING_LIST="missing_viral_only_samples.txt"
THREADS=12

mkdir -p "$HUMAN_ALIGN_DIR" "$VIRAL_ONLY_DIR"

# Read missing list into array
mapfile -t missing_combos < "$MISSING_LIST"

# Validate task ID
if [ -z "$SGE_TASK_ID" ] || [ "$SGE_TASK_ID" -gt "${#missing_combos[@]}" ]; then
    echo "SGE_TASK_ID out of bounds. Exiting."
    exit 1
fi

# Get sample-virus from list
base="${missing_combos[$((SGE_TASK_ID - 1))]}"
sample=${base%%__*}
virus=${base##*__}
viral_dir="$VIRAL_ALIGN_DIR/$base"

bam="$viral_dir/Aligned.sortedByCoord.out.bam"
viral_log="$viral_dir/Log.final.out"

if [[ ! -f "$bam" || ! -f "$viral_log" ]]; then
    echo "Missing BAM or STAR log for $base"
    exit 1
fi

# --- Viral STAR stats ---
unique_viral=$(grep "Uniquely mapped reads number" "$viral_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
multi_viral=$(grep "Number of reads mapped to multiple loci" "$viral_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
total_viral=$((unique_viral + multi_viral))

# --- Extract mapped viral reads ---
human_out="$HUMAN_ALIGN_DIR/$base"
mkdir -p "$human_out"

samtools view -b -F 4 "$bam" > "$human_out/mapped_viral.bam"
samtools sort -n "$human_out/mapped_viral.bam" -o "$human_out/mapped_viral.query.bam"
samtools fastq -1 "$human_out/mapped_viral_R1.fastq" -2 "$human_out/mapped_viral_R2.fastq" -0 /dev/null -s /dev/null -n "$human_out/mapped_viral.query.bam"

# --- Align to human genome ---
STAR --runThreadN $THREADS \
     --genomeDir "$HUMAN_INDEX" \
     --readFilesIn "$human_out/mapped_viral_R1.fastq" "$human_out/mapped_viral_R2.fastq" \
     --outFileNamePrefix "$human_out/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 9999

human_log="$human_out/Log.final.out"
human_bam="$human_out/Aligned.sortedByCoord.out.bam"

if [[ ! -f "$human_log" || ! -f "$human_bam" ]]; then
    echo "Human alignment failed for $base"
    exit 1
fi

# --- Human STAR stats ---
unique_human=$(grep "Uniquely mapped reads number" "$human_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
multi_human=$(grep "Number of reads mapped to multiple loci" "$human_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
total_human=$((unique_human + multi_human))

# --- Get human read names (with @ prefix) ---
viral_only_out="$VIRAL_ONLY_DIR/$base"
mkdir -p "$viral_only_out"
samtools view -F 4 "$human_bam" | cut -f1 | sort -u | sed 's/^/@/' > "$viral_only_out/human_mapped_readnames.txt"

# --- Properly filter viral FASTQ files using awk to skip 4-line blocks ---
filter_fastq() {
  input_fastq="$1"
  output_fastq="$2"
  grep_list="$3"

  awk -v list="$grep_list" '
    BEGIN {
      while ((getline line < list) > 0) {
        ids[line] = 1
      }
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

# --- Re-align viral-only reads to viral genome ---
realigned_out="$viral_only_out/realigned"
mkdir -p "$realigned_out"

STAR --runThreadN $THREADS \
     --genomeDir "${VIRAL_INDEX_BASE}/${virus}" \
     --readFilesIn "$viral_only_out/viral_only_R1.fastq" "$viral_only_out/viral_only_R2.fastq" \
     --outFileNamePrefix "$realigned_out/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMismatchNoverReadLmax 0.1 \
     --outFilterMismatchNmax 20

realigned_log="$realigned_out/Log.final.out"

# --- Stats from realignment ---
if [[ -f "$realigned_log" ]]; then
    viral_only_unique=$(grep "Uniquely mapped reads number" "$realigned_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
    viral_only_multi=$(grep "Number of reads mapped to multiple loci" "$realigned_log" | awk -F '|' '{gsub(/ /,"",$2); print $2}')
    viral_only_total=$((viral_only_unique + viral_only_multi))
else
    viral_only_unique=0
    viral_only_multi=0
    viral_only_total=0
fi

# --- Save stats ---
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

echo "Completed $base"