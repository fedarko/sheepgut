#! /usr/bin/env bash

OUTDIR=/Poppy/mfedarko/sheepgut/main-workflow/output
REF_FILE=$OUTDIR/all_edges.fasta
# This contains the "smaller" dataset; we also include the "larger" dataset
# (~5x as much data) in the alignment operation below.
SMALLER_READS_FILE=/Poppy/mkolmogo/sheep_meta/data/sheep_poop_CCS_dedup.fastq.gz

echo "Aligning reads to edges..."

# Use -x asm20 since these reads were generated using PacBio CCS; see
# https://github.com/lh3/minimap2#getting-started

# Use of multiple reads files based on comment from Heng Li here:
# https://github.com/lh3/minimap2/issues/191#issuecomment-399759935
minimap2 -ax asm20 $REF_FILE $SMALLER_READS_FILE /Poppy/mkolmogo/sheep_meta/data/ccs_sequel_II/*.fasta.gz > $OUTDIR/alignment.sam

echo "Did alignment."
