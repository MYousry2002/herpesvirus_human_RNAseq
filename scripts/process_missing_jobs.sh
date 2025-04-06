#!/bin/bash
#$ -N Align_herpesvirus_retry
#$ -cwd
#$ -o ../logs/retry_align_herpesvirus_$TASK_ID.out
#$ -e ../logs/retry_align_herpesvirus_$TASK_ID.err
#$ -pe smp 12
#$ -l h_vmem=64G
#$ -l h_rt=08:00:00
#$ -t 1-50

# Activate env
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Config
FASTQ_DIR="../data/GSE221091"
INDEX_DIR="../genomes/herpesviruses/STAR_indices"
OUTDIR="../workdir/alignments_herpesviruses"
THREADS=12
MISSING_LIST="./missing_jobs.txt"

# Get this task’s sample__virus
combo=$(sed -n "${SGE_TASK_ID}p" "$MISSING_LIST")
sample=${combo%%__*}
virus=${combo##*__}

fq1="${FASTQ_DIR}/${sample}_1.fastq"
fq2="${FASTQ_DIR}/${sample}_2.fastq"
index="${INDEX_DIR}/${virus}"
sample_out="${OUTDIR}/${sample}__${virus}"

# Skip if already done
if [[ -f "${sample_out}/Aligned.sortedByCoord.out.bam" ]]; then
    echo "Already done: $combo"
    exit 0
fi

# Run STAR
mkdir -p "$sample_out"
echo "🔬 Aligning $sample to $virus ..."

STAR --runThreadN $THREADS \
     --genomeDir "$index" \
     --readFilesIn "$fq1" "$fq2" \
     --outFileNamePrefix "${sample_out}/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFilterMismatchNoverReadLmax 0.1 \
     --outFilterMismatchNmax 20

echo "Finished $combo"