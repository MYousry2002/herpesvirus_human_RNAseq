#!/bin/bash

VIRUS_GENOMES="../genomes/herpesviruses"
INDEX_DIR="../genomes/herpesviruses/STAR_indices"
THREADS=12
SJDB_OVERHANG=84  # adjust based on your read length (e.g. read length - 1)

mkdir -p "$INDEX_DIR"

for fasta in "$VIRUS_GENOMES"/*.fasta; do
    virus_name=$(basename "$fasta" .fasta)
    output_dir="${INDEX_DIR}/${virus_name}"
    mkdir -p "$output_dir"
    
    echo "Building STAR index for $virus_name..."
    
    gtf="${VIRUS_GENOMES}/${virus_name}.gtf"
    
    if [ -f "$gtf" ]; then
        echo "Annotation file found for $virus_name: $gtf"
        STAR --runThreadN $THREADS \
             --runMode genomeGenerate \
             --genomeDir "$output_dir" \
             --genomeFastaFiles "$fasta" \
             --sjdbGTFfile "$gtf" \
             --sjdbGTFfeatureExon CDS \
             --sjdbOverhang $SJDB_OVERHANG \
             --genomeSAindexNbases 5
    else
        echo "No annotation file found for $virus_name, building index without annotation."
        STAR --runThreadN $THREADS \
             --runMode genomeGenerate \
             --genomeDir "$output_dir" \
             --genomeFastaFiles "$fasta" \
             --genomeSAindexNbases 5
    fi
done

echo "All viral STAR indices built at: $INDEX_DIR"