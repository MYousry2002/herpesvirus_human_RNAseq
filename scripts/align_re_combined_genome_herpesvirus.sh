#!/bin/bash
#$ -N Align_combined_herpesvirus
#$ -cwd
#$ -o ../logs/star_align_combined_herpesvirus_$TASK_ID.out
#$ -e ../logs/star_align_combined_herpesvirus_$TASK_ID.err
#$ -pe smp 8
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
#$ -t 1-250

# If SGE_TASK_ID is not set (e.g. when running interactively), set it to 1
if [ -z "$SGE_TASK_ID" ]; then
    export SGE_TASK_ID=1
fi

# Activate environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Paths
FASTQ_DIR="../data/GSE221091"
INDEX_DIR="../genomes/herpesviruses/STAR_indices/combined_herpes_viruses"
OUTDIR="../workdir/alignment_herpesvirus"
THREADS=12
SAMPLE_LIST="../data/GSE221091/sample_list.txt"

# Get sample name
sample=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')
echo "Task ${SGE_TASK_ID} processing sample: [$sample]"

if [ -z "$sample" ]; then
    echo "No sample found for task ${SGE_TASK_ID}. Exiting."
    exit 1
fi

fq1="$FASTQ_DIR/${sample}_1.fastq"
fq2="$FASTQ_DIR/${sample}_2.fastq"

if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
    echo "Missing FASTQ files for $sample"
    exit 1
fi

# Define combined sample output
virus="combined_herpes_viruses"
sample_out="${OUTDIR}/${sample}__${virus}"

# Clean existing results if present
if [[ -d "$sample_out" ]]; then
    echo "Removing previous results for $sample with $virus"
    rm -rf "$sample_out"
fi

mkdir -p "$sample_out"

echo "Aligning $sample to $virus..."

STAR --runThreadN $THREADS \
     --genomeDir "$INDEX_DIR" \
     --readFilesIn "$fq1" "$fq2" \
     --outFileNamePrefix "${sample_out}/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMismatchNoverReadLmax 0.1 \
     --outFilterMismatchNmax 20 \
     --outFilterMultimapNmax 9999


echo "Finished re-aligning $sample to $virus"