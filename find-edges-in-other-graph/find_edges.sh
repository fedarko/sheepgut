#! /usr/bin/env bash

# File of the edge sequences in the smaller graph to search for in the
# larger graph
REF_FILE=/Poppy/mfedarko/sheep_metagenome/166_6018_7998.fasta
# Sequences in the larger graph (we're using the non-haplotype-mode graph,
# although the smaller graph did use the haplotype mode)
READS_FILE=/Poppy/mkolmogo/sheep_meta/flye_big_2.8/assembly.fasta

# Use -x asm20 since these reads were generated using PacBio CCS; see
# https://github.com/lh3/minimap2#getting-started
minimap2 -ax asm20 $REF_FILE $READS_FILE > graph_edges_to_166_etal_nonhaplo.sam
