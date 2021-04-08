#! /usr/bin/env bash

REF_FILE=newseqs.fasta
ALN_FILE=aln-sorted.bam
OUT_FILE=pileup.txt

# We want to compute the pileups for just a few edges.
samtools mpileup $ALN_FILE -f $REF_FILE > $OUT_FILE
