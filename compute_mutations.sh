#! /usr/bin/env bash

REF_FILE=/Poppy/mfedarko/sheep_metagenome/selected_edges.fasta
ALN_FILE=/Poppy/mfedarko/sheep_metagenome/saln-scf-sorted.bam
OUT_FILE=/Poppy/mfedarko/sheep_metagenome/scf-pileup.txt

# We want to compute the pileups for just a few edges.
samtools mpileup $ALN_FILE -f $REF_FILE > $OUT_FILE
