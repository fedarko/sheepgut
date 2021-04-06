#! /usr/bin/env bash

REF_FILE=/Poppy/mfedarko/sheep_metagenome/redo_work_big_nonhaplo_graph/all_edges.fasta
READS_FILE=/Poppy/mkolmogo/sheep_meta/data/sheep_poop_CCS_dedup.fastq.gz

# Use -x asm20 since these reads were generated using PacBio CCS; see
# https://github.com/lh3/minimap2#getting-started

# Use of multiple reads files based on comment from Heng Li here:
# https://github.com/lh3/minimap2/issues/191#issuecomment-399759935
minimap2 -ax asm20 $REF_FILE $READS_FILE /Poppy/mkolmogo/sheep_meta/data/ccs_sequel_II/*.fasta.gz > alignment.sam
