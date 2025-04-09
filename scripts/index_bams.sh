#!/bin/bash

for bam in ../workdir/viral_only/*/realigned/Aligned.sortedByCoord.out.bam; do
    samtools index "$bam"
done