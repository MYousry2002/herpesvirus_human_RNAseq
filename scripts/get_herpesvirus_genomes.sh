#!/bin/bash

# Set target directory
TARGET_DIR="../genomes/herpesviruses"

# Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Change into the target directory
cd "$TARGET_DIR" || exit

# Download each herpesvirus genome as FASTA
efetch -db nucleotide -format fasta -id KT899744.1 > HSV1_KOS.fasta              # Herpes Simplex Virus 1
efetch -db nucleotide -format fasta -id OM370995.1 > HSV2_StrainG.fasta         # Herpes Simplex Virus 2
efetch -db nucleotide -format fasta -id KU926311.1 > VZV_Ellen.fasta            # Varicella Zoster Virus
efetch -db nucleotide -format fasta -id V01555.2 > EBV.fasta                    # Epstein-Barr Virus
efetch -db nucleotide -format fasta -id FJ616285.1 > HCMV.fasta                 # Human Cytomegalovirus
efetch -db nucleotide -format fasta -id MW536483.1 > HHV6B.fasta                # Human Herpesvirus 6B
efetch -db nucleotide -format fasta -id MZ712172.1 > KSHV.fasta                 # Kaposi’s Sarcoma Virus (HHV-8)
efetch -db nucleotide -format fasta -id AF037218.1 > HHV7.fasta                 # Human Herpesvirus 7

# Combine them all into one file
cat *.fasta > combined_herpes_viruses.fasta

echo "All herpesvirus genomes downloaded and combined into:"
echo "$(realpath combined_herpes_viruses.fasta)"