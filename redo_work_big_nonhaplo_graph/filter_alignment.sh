#! /usr/bin/env bash

# try to remove "peaks" in coverage plots (alignment artifacts)
# by just retaining primary alignments.
# based on https://github.com/lh3/minimap2/issues/416#issuecomment-499095442.
# We currently ONLY filter out secondary alignments, but not supplementary
# alignments; either possible.
samtools view -b -F0x100 alignment.sam > alignment.bam
