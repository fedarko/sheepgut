#! /usr/bin/env bash
#
# To see what's up with the deletions in BACT2 and coverage drops in CAMP +
# BACT1, let's rerun alignment but keep all BAM files (including before OSA
# and partially-mapped-read filtering). This should tell us either that 1)
# these filtering steps *cause* these patterns, or 2) they don't. (Or something
# in the middle I guess)
#
# The filename of this script starts with a "5" but it isn't really necessary
# to run the other stuff before this lol (although you should've run the main
# workflow stuff, b/c we rely on all_edges.fasta existing... although at that
# point you could just use gfa-to-fasta to make it so whatevs)

CONTIGS=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta
GFA=/Poppy/mfedarko/misc-data/sheepgut_flye_big_2.8_graph.gfa
read READS < ../../config/input-reads

# (the --no-rm-tmp-bam thing is the important part!)
strainFlye align \
    --contigs $CONTIGS \
    --graph $GFA \
    --no-rm-tmp-bam \
    --verbose \
    --output-dir output/aln \
    $READS
