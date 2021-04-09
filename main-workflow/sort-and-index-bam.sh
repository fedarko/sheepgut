#! /usr/bin/env bash

# This program takes two command-line arguments: an input BAM file, and
# the output path to write the sorted BAM file to (and for which an index
# should be created).

echo "Sorting and indexing BAM file $1..."

# https://stackoverflow.com/a/638980
if [ ! -f $1 ]; then
    echo "File $1 not found?"
    # Return a nonzero value to indicate an error: from
    # https://stackoverflow.com/a/50265513
    exit 1
fi

# Based on http://quinlanlab.org/tutorials/samtools/samtools.html
samtools sort $1 -T output/_tmpfile -o $2

echo "Sorted the BAM file $1: the sorted BAM file is named $2."

echo "Indexing BAM file $2..."

samtools index $2

echo "Indexed the sorted BAM file $2."
