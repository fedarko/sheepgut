#! /usr/bin/env bash

# Avoid samtools sort creating hundreds of temporary files in the CWD
# UPDATE: never mind this samtools version (0.1.19...) is old enough that it
# doesn't support -T, so we'll just live with tmp files being in this directory
# TMPDIR=/Poppy/mfedarko/sheep_metagenome/tmp_samtools_files
#SAMFILE=alignment-filtered.sam
BAMFILE=alignment.bam

#samtools view -bS $SAMFILE > $BAMFILE
#
#echo "did view"
# # Based on http://quinlanlab.org/tutorials/samtools/samtools.html
samtools sort $BAMFILE -T _tmpfile -o aln-sorted.bam

echo "did sort"

samtools index aln-sorted.bam

echo "did index"
