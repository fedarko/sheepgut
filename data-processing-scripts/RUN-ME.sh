#! /usr/bin/env bash

# NOTE: Since the alignment process works against ALL edges, the only part
# that needs to be rerun to examine different edge sequences is the stuff
# starting at ./filter-partially-mapped-reads.py (and stuff in
# bam2pointinfojsons.py; I should ideally make the edge sequences configurable
# from in this script or something)

# NOTE 2: The mut-analyses conda env should be activated before running
# this, since having later version of samtools (1.7 as of writing) plus
# pysam/etc. is needed

./gfa-to-fasta.py
echo "Created FASTA of all edges in the assembly graph."

echo "Aligning reads to edges..."
./align-reads-to-edges.sh
echo "Did alignment."

./filter-secondary-alignments.sh
echo "Filtered secondary alignments and created a BAM file."

./sort-and-index-bam.sh alignment.bam aln-sorted.bam
echo "Sorted and indexed the BAM file."

./filter-partially-mapped-reads.py
echo "Filtered to fully (ish) aligned reads."

# I'm not sure if sorting this BAM again is needed, but it shouldn't cause
# problems to re-sort it anyway (beyond being somewhat inefficient).
# We do need to index it so that we can use pileup in bam2pointinfojsons.py.
./sort-and-index-bam.sh pmread-filtered-aln.bam fully-filtered-and-sorted-aln.bam

./bam-to-jsons.py
echo "Done."
