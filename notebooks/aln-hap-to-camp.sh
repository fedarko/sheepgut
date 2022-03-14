#!/usr/bin/env bash

CAMP=/Poppy/mfedarko/sheepgut/seqs/edge_6104.fasta
PREFIX=phasing-data/smoothed-reads/camp_g1217
OUTDIR=$PREFIX/aln

# Yeah this could be done in a simple for loop but bash scares me

# metaFlye
minimap2 -ax asm20 --MD $CAMP $PREFIX/metaflye/assembly.fasta > $OUTDIR/metaflye.sam
minimap2 -ax asm20 --MD $CAMP \
    $PREFIX/metaflye_keephaps/assembly.fasta > $OUTDIR/metaflye_keephaps.sam
    
# jumboDBG
minimap2 -ax asm20 --MD $CAMP \
    $PREFIX/jumbodbg_k5001/graph.fasta > $OUTDIR/jumbodbg.sam
    
# LJA
minimap2 -ax asm20 --MD $CAMP \
    $PREFIX/lja_withec/assembly.fasta > $OUTDIR/lja_withec.sam
minimap2 -ax asm20 --MD $CAMP \
    $PREFIX/lja_noec/assembly.fasta > $OUTDIR/lja_noec.sam
minimap2 -ax asm20 --MD $CAMP \
    $PREFIX/lja_covfilt_10x/assembly.fasta > $OUTDIR/lja_covfilt_10x.sam
