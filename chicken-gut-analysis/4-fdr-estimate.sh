#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

# Let's use all currently-available decoy contexts, just for the sake of
# demonstration. (will eventually add an "all of the above" option...)
strainFlye fdr estimate \
    --contigs $CGDIR/sf/contigs.fasta \
    --bam $CGDIR/sf/aln/final.bam \
    --bcf $CGDIR/sf/call-r-minr3/naive-calls.bcf \
    --diversity-indices $CGDIR/sf/call-r-minr3/diversity-indices.tsv \
    --decoy-contexts Full \
    --decoy-contexts CP2 \
    --decoy-contexts Tv \
    --decoy-contexts Nonsyn \
    --decoy-contexts Nonsense \
    --decoy-contexts CP2Tv \
    --decoy-contexts CP2Nonsyn \
    --decoy-contexts CP2Nonsense \
    --decoy-contexts TvNonsyn \
    --decoy-contexts TvNonsense \
    --decoy-contexts CP2TVNonsense \
    --decoy-min-average-coverage 100 \
    --output-dir $CGDIR/sf/fdr-estimate-dmac100x/
