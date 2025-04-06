#!/bin/bash
#$ -N Align_herpesvirus
#$ -cwd
#$ -o ../logs/star_align_herpesvirus_$TASK_ID.out
#$ -e ../logs/star_align_herpesvirus_$TASK_ID.err
#$ -pe smp 12
#$ -l h_vmem=64G
#$ -l h_rt=24:00:00
#$ -t 1-250

# If SGE_TASK_ID is not set (e.g. when running interactively), set it to 1
if [ -z "$SGE_TASK_ID" ]; then
    export SGE_TASK_ID=1
fi

# Activate environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Define paths using absolute paths for reliability
FASTQ_DIR="../data/GSE221091"
INDEX_DIR="../genomes/herpesviruses/STAR_indices"
OUTDIR="../workdir/alignment_herpesvirus"
THREADS=12
SAMPLE_LIST="../data/GSE221091/sample_list.txt"

# Get sample for this task
sample=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')
echo "Task ${SGE_TASK_ID} processing sample: [$sample]"

# Check if sample variable is empty
if [ -z "$sample" ]; then
    echo "No sample found for task ${SGE_TASK_ID}. Exiting."
    exit 1
fi

fq1="$FASTQ_DIR/${sample}_1.fastq"
fq2="$FASTQ_DIR/${sample}_2.fastq"

# Check if both FASTQ files exist; if missing, skip this sample.
if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
    echo "Missing FASTQ files for $sample"
    exit 1
fi

# Loop over all viral genome STAR indices
for index in "$INDEX_DIR"/*; do
    virus=$(basename "$index")
    sample_out="${OUTDIR}/${sample}__${virus}"
    mkdir -p "$sample_out"
    
    # Check if this sample-viral index pair is already processed
    if [[ -f "${sample_out}/Aligned.sortedByCoord.out.bam" ]]; then
        echo "Already processed $sample for virus $virus. Skipping..."
        continue
    fi

    echo "Aligning $sample to $virus ..."
    
    STAR --runThreadN $THREADS \
         --genomeDir "$index" \
         --readFilesIn "$fq1" "$fq2" \
         --outFileNamePrefix "${sample_out}/" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outFilterMismatchNoverReadLmax 0.1 \
         --outFilterMismatchNmax 20
done

echo "Finished processing sample: $sample"