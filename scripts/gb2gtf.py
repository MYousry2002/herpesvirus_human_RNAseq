#!/usr/bin/env python
"""
A simple GenBank-to-GTF converter using Biopython.
Usage: python gb2gtf.py input.gb output.gtf

This script produces a basic GTF file by iterating through the features
in the GenBank file. It outputs the feature type, start, end, strand, and
an attribute field containing gene information (if available).
"""

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    sys.exit("Usage: python gb2gtf.py input.gb output.gtf")

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as out:
    out.write("##gtf-version 2.2\n")
    for record in SeqIO.parse(input_file, "genbank"):
        seqid = record.id
        for feature in record.features:
            if feature.type == "source":
                continue
            start = int(feature.location.start) + 1  # GTF is 1-based
            end = int(feature.location.end)
            strand = "+"
            if feature.location.strand == -1:
                strand = "-"
            elif feature.location.strand == 0 or feature.location.strand is None:
                strand = "."
            score = "."
            frame = "."
            # Build an attribute string: if gene or locus_tag is present.
            if "gene" in feature.qualifiers:
                gene = feature.qualifiers["gene"][0]
                attr = f'gene_id "{gene}"; gene_name "{gene}";'
            elif "locus_tag" in feature.qualifiers:
                tag = feature.qualifiers["locus_tag"][0]
                attr = f'gene_id "{tag}"; gene_name "{tag}";'
            else:
                attr = f'gene_id "{feature.type}";'
            out.write(f"{seqid}\tBiopython\t{feature.type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr}\n")