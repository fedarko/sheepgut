#! /usr/bin/env bash

# Just runs LoFreq on the three selected MAGs (CAMP, BACT1, BACT2)
lofreq call \
    -f output/selected-mags.fasta \
    -o output/lofreq.vcf \
    output/selected-mags.bam
