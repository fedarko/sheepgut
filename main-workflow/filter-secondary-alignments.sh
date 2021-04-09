#! /usr/bin/env bash

OUTDIR=/Poppy/mfedarko/sheepgut/data-processing-scripts/output

echo "Filtering secondary alignments and converting to BAM..."

# Try to remove "peaks" in coverage plots (alignment artifacts)
# by just retaining primary alignments.
#
# Based on https://github.com/lh3/minimap2/issues/416#issuecomment-499095442.
#
# We currently just filter out secondary alignments, but not supplementary
# alignments; we leave those in so we can use those while making decisions
# about partially-mapped read filtering in a later step.
samtools view -b -F0x100 $OUTPUT/alignment.sam > $OUTPUT/alignment.bam

echo "Done."
