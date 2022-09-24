#! /usr/bin/env bash

CONTIGS=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta
BAM=/Poppy/mfedarko/sheepgut/main-workflow/output/fully-filtered-and-sorted-aln.bam

strainFlye fdr estimate \
    --contigs $CONTIGS \
    --bam $BAM \
    --bcf output/call-r3/naive-calls.bcf \
    --decoy-contig "edge_6104" \
    --decoy-contexts "Everything" \
    --output-dir output/fdr-estimate-r
