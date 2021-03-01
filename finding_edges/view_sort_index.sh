#! /usr/bin/env bash

PREFIX=/Poppy/mfedarko/sheep_metagenome/finding_edges
SAMFILE=$PREFIX/graph_edges_to_166_etal.sam
BAMFILE=$PREFIX/graph_edges_to_166_etal.bam

samtools view -bS $SAMFILE > $BAMFILE

echo "did view"
# # Based on http://quinlanlab.org/tutorials/samtools/samtools.html
# # The -o option didn't seem to work -- samtools just wrote a bunch of junk
# # to stdout regardless -- so I'm just manually piping the output into a new
# # BAM file.
samtools sort $BAMFILE aln-sorted > $PREFIX/alignment_sorted.bam

echo "did sort"

samtools index $PREFIX/aln-sorted.bam

echo "did index"
