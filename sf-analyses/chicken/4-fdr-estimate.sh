#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

# Let's use all currently-available decoy contexts, just for the sake of
# demonstration. (will eventually add an "all of the above" option...)
strainFlye fdr estimate \
    --contigs $CGDIR/sf/contigs.fasta \
    --bam $CGDIR/sf/aln/final.bam \
    --bcf $CGDIR/sf/call-p-minp100/naive-calls.bcf \
    --diversity-indices $CGDIR/sf/call-p-minp100/diversity-indices.tsv \
    --decoy-contexts Everything \
    --decoy-min-average-coverage 100 \
    --output-dir $CGDIR/sf/fdr-estimate-dmac100x-p/
