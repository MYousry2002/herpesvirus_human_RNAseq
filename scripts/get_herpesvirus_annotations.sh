#!/bin/bash
# This script downloads GenBank records for a set of viruses,
# converts them to GTF format using a custom Python script (gb2gtf.py),
# and stores the GTF files in ../genomes/herpesviruses.

# Create output directory if it doesn't exist
OUTDIR="../genomes/herpesviruses"
mkdir -p "$OUTDIR"

# Define an array with entries in the format:
# "VirusName;Accession;BaseFilename"
# Here the third field is the desired (short) name for the output file.
entries=(
  "Herpes_Simplex_1_KOS;KT899744.1;HSV1_KOS"
  "Herpes_Simplex_2_Strain_G;OM370995.1;HSV2_StrainG"
  "Varicella_Zoster_Virus_Ellen;KU926311.1;VZV_Ellen"
  "Epstein_Barr_Virus;V01555.2;EBV"
  "Human_cytomegalovirus;FJ616285.1;HCMV"
  "HHV6B;MW536483.1;HHV6B"
  "Kaposi_Sarcoma_HHV8;MZ712172.1;KSHV"
  "Human_Herpes_7;AF037218.1;HHV7"
)

# Loop over each entry and process the virus
for entry in "${entries[@]}"; do
    IFS=';' read -r virus accession basefilename <<< "$entry"
    echo "Downloading GenBank record for $virus (Accession: $accession)..."
    efetch -db nucleotide -format gb -id "$accession" > "${basefilename}.gb"
    if [ $? -ne 0 ]; then
      echo "Error downloading $accession"
      continue
    fi
    echo "Downloaded ${basefilename}.gb"
    
    # Convert the GenBank file to GTF using the custom Python script
    python gb2gtf.py "${basefilename}.gb" "${basefilename}.gtf"
    if [ $? -ne 0 ]; then
      echo "Error converting ${basefilename}.gb to GTF"
      continue
    fi
    echo "Converted to ${basefilename}.gtf"
    
    # Move the final GTF file to the output directory
    mv "${basefilename}.gtf" "$OUTDIR"
    
    # Optionally, remove the intermediate GenBank file
    rm -f "${basefilename}.gb"
    
    echo "Annotation for $virus saved to $OUTDIR/${basefilename}.gtf"
done

echo "All virus annotations processed."