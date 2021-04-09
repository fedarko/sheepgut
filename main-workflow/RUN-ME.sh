#! /usr/bin/env bash

# NOTE: Since the alignment process works against ALL edges, if it's desired
# to rerun this analysis while focusing on different edges in the graph,
# the only parts that need to be rerun are the stuff starting at and below
# ./filter-partially-mapped-reads.py.
#
# Ideally I'd make the edge sequences configurable in this script or something.
#
# NOTE 2: The mut-analyses conda env should be activated before running
# this, since having later version of samtools (1.7 as of writing) plus
# pysam/etc. is needed

OUTDIR=/Poppy/mfedarko/sheepgut/main-workflow/output

echo "Starting main analysis workflow."

./gfa-to-fasta.py

./align-reads-to-edges.sh

./filter-secondary-alignments.sh

# We need to have indexed the BAM file in order to use fetch() in Pysam
# in the partially-mapped reads filtering script below. (And, at least
# as far as I'm aware, we need to sort a BAM file before indexing it.)
./sort-and-index-bam.sh $OUTDIR/alignment.bam $OUTDIR/aln-sorted.bam

./filter-partially-mapped-reads.py

# I'm not sure if sorting this BAM again is strictly necessary, but we err
# on the side of safety and re-sort it anyway. We do need to index it so
# that we can use pileup in the BAM-to-JSON script.
./sort-and-index-bam.sh $OUTDIR/pmread-filtered-aln.bam $OUTDIR/fully-filtered-and-sorted-aln.bam

./bam-to-jsons.py

echo "Done with the main analysis workflow!"
