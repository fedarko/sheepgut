#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

# We'll use the sequences from the .p_ctg.gfa file; these correspond to the
# "final" assembled sequences from hifiasm-meta, according to
# https://github.com/xfengnefx/hifiasm-meta/issues/10.

strainFlye utils gfa-to-fasta \
    --graph $CGDIR/asm/chicken.hifiasm-meta.p_ctg.gfa \
    --output-fasta $CGDIR/sf/contigs.fasta

# Align reads to these hifiasm-meta contigs

strainFlye align \
    --contigs $CGDIR/sf/contigs.fasta \
    --graph $CGDIR/asm/chicken.hifiasm-meta.p_ctg.gfa \
    --output-dir $CGDIR/sf/aln/
