#! /usr/bin/env bash

CONTIGS=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta
BAM=/Poppy/mfedarko/sheepgut/main-workflow/output/fully-filtered-and-sorted-aln.bam

# Call @ r = 3
strainFlye call r-mutation \
    --contigs $CONTIGS \
    --bam $BAM \
    --min-r 3 \
    --div-index-r-list "2,3,4,5,10,20,50,100,250,500" \
    --verbose \
    --output-dir output/call-r3
