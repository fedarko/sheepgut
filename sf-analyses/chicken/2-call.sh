#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

# TODO in the future, would be a good idea to lower --min-read-number
# to allow for easier target contig selection for the FDR curve plots
# (and also, probably better decoy contig selection)

strainFlye call p-mutation \
    --contigs $CGDIR/sf/contigs.fasta \
    --bam $CGDIR/sf/aln/final.bam \
    --min-p 100 \
    --verbose \
    --output-dir $CGDIR/sf/call-p-minp100/
