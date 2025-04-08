# In a working script dir:
SAMPLE_LIST=../data/GSE221091/sample_list.txt
INDEX_DIR=../genomes/herpesviruses/STAR_indices
OUTDIR=../workdir/alignments_herpesviruses

# Get virus names
ls $INDEX_DIR | sort > virus_list.txt

# Create list of missing sample__virus combinations
> missing_jobs.txt  # reset the file

while read sample; do
    while read virus; do
        outfile="${OUTDIR}/${sample}__${virus}/Aligned.sortedByCoord.out.bam"
        if [[ ! -f "$outfile" ]]; then
            echo "${sample}__${virus}" >> missing_jobs.txt
        fi
    done < virus_list.txt
done < $SAMPLE_LIST