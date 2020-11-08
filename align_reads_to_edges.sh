#! /usr/bin/env bash

REF_FILE=/Poppy/mfedarko/sheep_metagenome/selected_edges.fasta
READS_FILE=/Poppy/mkolmogo/sheep_meta/data/sheep_poop_CCS_dedup.fastq.gz

# Use -x asm20 since these reads were generated using PacBio CCS; see
# https://github.com/lh3/minimap2#getting-started
minimap2 -ax asm20 $REF_FILE $READS_FILE /Poppy/mkolmogo/sheep_metadata/data/ccs_sequel_II/*.fasta.gz > 5xalignment.sam
