#! /usr/bin/env bash

# Use of BED file to filter to the three selected MAGs
# ("references"/"chromosomes", going by samtools parlance) from
# http://seqanswers.com/forums/showthread.php?t=49476
samtools view -b -L selected-mags.bed output/fully-filtered-and-sorted-aln.bam > selected-mags.bam
