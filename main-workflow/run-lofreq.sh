#! /usr/bin/env bash

lofreq call \
    -f output/all_edges.fasta \
    -o output/lofreq.vcf \
    output/fully-filtered-and-sorted-aln.bam
