#! /usr/bin/env bash

# Avoid samtools sort creating hundreds of temporary files in the CWD
# UPDATE: never mind this samtools version (0.1.19...) is old enough that it
# doesn't support -T, so we'll just live with tmp files being in this directory
# TMPDIR=/Poppy/mfedarko/sheep_metagenome/tmp_samtools_files
PREFIX=/Poppy/mfedarko/sheep_metagenome

samtools view -bS $PREFIX/selected_alignment_less_sc_filtered.sam > $PREFIX/selected_alignment_less_sc_filtered.bam

echo "did view"
# # Based on http://quinlanlab.org/tutorials/samtools/samtools.html
# # The -o option didn't seem to work -- samtools just wrote a bunch of junk
# # to stdout regardless -- so I'm just manually piping the output into a new
# # BAM file.
samtools sort $PREFIX/selected_alignment_less_sc_filtered.bam saln-scf-sorted > $PREFIX/selected_alignment_scf_sorted.bam

echo "did sort"

samtools index $PREFIX/saln-scf-sorted.bam

echo "did index"
