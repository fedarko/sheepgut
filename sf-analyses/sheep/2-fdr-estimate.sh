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

strainFlye fdr estimate \
    --contigs $CONTIGS \
    --bam $BAM \
    --bcf output/call-p15/naive-calls.bcf \
    --decoy-contig "edge_6104" \
    --decoy-contexts "Everything" \
    --output-dir output/fdr-estimate
