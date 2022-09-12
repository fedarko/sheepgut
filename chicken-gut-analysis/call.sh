#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

strainFlye call r-mutation \
    --contigs $CGDIR/sf/contigs.fasta \
    --bam $CGDIR/sf/aln/final.bam \
    --min-r 3 \
    --div-index-r-list "2,3,4,5,10,20,50,100" \
    --verbose \
    --output-dir $CGDIR/sf/call-r-minr3/
