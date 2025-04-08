#!/bin/bash

GENOME_DIR="../genomes/herpesviruses"
TMP_DIR="$GENOME_DIR/tmp_combined"
COMBINED_FASTA="$GENOME_DIR/combined_herpes_viruses.fasta"
COMBINED_GTF="$GENOME_DIR/combined_herpes_viruses.gtf"
STAR_INDEX_DIR="$GENOME_DIR/STAR_indices/combined_herpes_viruses"
THREADS=12

mkdir -p "$TMP_DIR"

echo "Cleaning and prefixing FASTA and GTF..."

# Clean and prefix each virus's files
for fasta in "$GENOME_DIR"/*.fasta; do
    virus=$(basename "$fasta" .fasta)

    # Fix FASTA headers
    awk -v prefix="${virus}__" '/^>/ {$0=">"prefix substr($0,2)} 1' "$fasta" \
        > "$TMP_DIR/${virus}.prefixed.fasta"
done

for gtf in "$GENOME_DIR"/*.gtf; do
    virus=$(basename "$gtf" .gtf)

    # Fix GTF contig names
    awk -v prefix="${virus}__" '{if($0 !~ /^#/){$1=prefix $1} print}' OFS="\t" "$gtf" \
        > "$TMP_DIR/${virus}.prefixed.gtf"
done

echo "Merging all cleaned FASTA files..."
cat "$TMP_DIR"/*.prefixed.fasta > "$COMBINED_FASTA"

echo "Merging all cleaned GTF files..."
cat "$TMP_DIR"/*.prefixed.gtf > "$COMBINED_GTF"

echo "Preview of FASTA headers:"
grep '^>' "$COMBINED_FASTA" | head

echo "Preview of GTF contigs:"
cut -f1 "$COMBINED_GTF" | sort | uniq | head

echo "Building STAR genome index..."
mkdir -p "$STAR_INDEX_DIR"
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$COMBINED_FASTA" \
     --sjdbGTFfile "$COMBINED_GTF" \
     --sjdbGTFfeatureExon CDS \
     --sjdbOverhang 85

echo "STAR index created in: $STAR_INDEX_DIR"

# clean up
rm -r "$TMP_DIR"