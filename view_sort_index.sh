#! /usr/bin/env bash

# Avoid samtools sort creating hundreds of temporary files in the CWD
# UPDATE: never mind this samtools version (0.1.19...) is old enough that it
# doesn't support -T, so we'll just live with tmp files being in this directory
# TMPDIR=/Poppy/mfedarko/sheep_metagenome/tmp_samtools_files
PREFIX=/Poppy/mfedarko/sheep_metagenome
SAMFILE=$PREFIX/5xalignment_less_soft_clipping.sam
BAMFILE=$PREFIX/5xalignment_less_soft_clipping.bam

samtools view -bS $SAMFILE > $BAMFILE

echo "did view"
# # Based on http://quinlanlab.org/tutorials/samtools/samtools.html
# # The -o option didn't seem to work -- samtools just wrote a bunch of junk
# # to stdout regardless -- so I'm just manually piping the output into a new
# # BAM file.
samtools sort $BAMFILE saln-sorted > $PREFIX/selected_alignment_sorted.bam

echo "did sort"

samtools index $PREFIX/saln-sorted.bam

echo "did index"
