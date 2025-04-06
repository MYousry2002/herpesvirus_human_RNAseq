#!/bin/bash

# Create and move into the genome directory
TARGET_DIR="../genomes/human"
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR" || exit

# Download latest GRCh38 (primary assembly, soft-masked) from Ensembl
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

# decompress the files
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

echo "Human genome and annotations downloaded to: $(pwd)"