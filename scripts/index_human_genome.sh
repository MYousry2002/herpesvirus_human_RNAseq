#!/bin/bash
#$ -N STARindex_human
#$ -cwd
#$ -o ../logs/star_index_human_$JOB_ID.out
#$ -e ../logs/star_index_human_$JOB_ID.err
#$ -pe smp 12
#$ -l h_vmem=4G
#$ -l h_rt=08:00:00
#$ -m e

# Load conda environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate herpesvirus

# Define paths
HUMAN_FASTA="../genomes/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
HUMAN_GTF="../genomes/human/Homo_sapiens.GRCh38.110.gtf"
HUMAN_INDEX="../genomes/human/STAR_index"

# Create output directory
mkdir -p "$HUMAN_INDEX"

# Build STAR index
STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir "$HUMAN_INDEX" \
     --genomeFastaFiles "$HUMAN_FASTA" \
     --sjdbGTFfile "$HUMAN_GTF" \
     --sjdbOverhang 100

echo "STAR index build finished."