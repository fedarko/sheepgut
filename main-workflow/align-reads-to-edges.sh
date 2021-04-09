#! /usr/bin/env bash

REF_FILE=output/all_edges.fasta

# Locate the input read files. Use of "read" (the bash command, not the
# sequencing data!) based on http://mywiki.wooledge.org/BashFAQ/028.
#
# This should include the full dataset (i.e. the smaller dataset and the
# additional ~5x-as-much sequencing data). Globs are OK to include in this
# variable, since this will likely be spread across many FASTA / FASTQ files;
# the globs will be evaluated when we expand this variable as $READS_FILES --
# see https://unix.stackexchange.com/a/314810.
read READS_FILES < ../config/input-reads
echo "The input reads files are $READS_FILES."
echo "(After expansion, it looks like there are a total of `echo $READS_FILES | wc -w` files included.)"

echo "Aligning reads to edges..."

# Use -x asm20 since these reads were generated using PacBio CCS; see
# https://github.com/lh3/minimap2#getting-started

# Use of multiple reads files based on comment from Heng Li here:
# https://github.com/lh3/minimap2/issues/191#issuecomment-399759935
minimap2 -ax asm20 $REF_FILE $READS_FILES > output/alignment.sam

echo "Did alignment."
