#! /usr/bin/env bash

CGDIR=/Poppy/mfedarko/chicken-gut-meta

# We'll use the sequences from the .p_ctg.gfa file; these correspond to the
# "final" assembled sequences from hifiasm-meta, according to
# https://github.com/xfengnefx/hifiasm-meta/issues/10.

strainFlye utils gfa-to-fasta \
    --graph $CGDIR/asm/chicken.hifiasm-meta.p_ctg.gfa \
    --output-fasta $CGDIR/sf/contigs.fasta

# Align reads to these hifiasm-meta contigs
# (note that the reads come at the end of the command, since we allow the user
# to specify a variable number of reads files -- although in this case, it's
# just one FASTQ file)

strainFlye align \
    --contigs $CGDIR/sf/contigs.fasta \
    --graph $CGDIR/asm/chicken.hifiasm-meta.p_ctg.gfa \
    --verbose \
    --output-dir $CGDIR/sf/aln/ \
    $CGDIR/reads/SRR15214153.fastq
