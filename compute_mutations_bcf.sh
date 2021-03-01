#! /usr/bin/env bash

REF_FILE=/Poppy/mfedarko/sheep_metagenome/selected_edges.fasta
ALN_FILE=/Poppy/mfedarko/sheep_metagenome/saln-sorted.bam
OUT_FILE=/Poppy/mfedarko/sheep_metagenome/5xdata-pileup-gI.bcf

# Pileups, but:
# -g output to BCF
# -I ignore indels (https://www.biostars.org/p/112255/)
samtools mpileup $ALN_FILE -f $REF_FILE -g -I > $OUT_FILE
