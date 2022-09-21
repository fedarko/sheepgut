#! /usr/bin/env bash
#
# Performs mutation calling using the actual strainFlye pipeline code (not the
# messy ad hoc code used in notebooks/, for example).
#
# This enables us to use some of the more recently developed strainFlye
# functionality on the SheepGut dataset. For example, drawing some recently
# implemented types of FDR curves, or computing a p-value for the longest gap
# between mutations in a contig.

CONTIGS=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta
BAM=/Poppy/mfedarko/sheepgut/main-workflow/output/fully-filtered-and-sorted-aln.bam

# Call @ p = 0.15% (used for FDR curve plotting)
# (We specify a custom set of values of p for the diversity index just so we
# can replicate the Diversity Index figure in the paper while we're at it)
strainFlye call p-mutation \
    --contigs $CONTIGS \
    --bam $BAM \
    --min-p 15 \
    --div-index-p-list "25,50,100,200,500,1000,2500,5000" \
    --verbose \
    --output-dir output/call-p15

# Call @ p = 0.50% (used for computing coldspots, to match the set of mutations
# we used to identify coldspots in the submitted version of the paper)
strainFlye call p-mutation \
    --contigs $CONTIGS \
    --bam $BAM \
    --min-p 50 \
    --verbose \
    --output-dir output/call-p50
