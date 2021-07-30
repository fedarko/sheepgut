#! /usr/bin/env bash
# Runs the "main workflow" scripts, detailed in more depth in the README.
#
# This takes a while (longest step will likely be the second, in which reads
# are aligned to all edges in the assembly graph).

# Attempt to create an output directory if it doesn't exist. The "output" name
# for this directory is used a lot throughout this repository, so I do not
# recommend changing this. (We could make it configurable, like the input file
# paths, but I'm not sure that would be worth the effort.)
OUTDIR=output/
mkdir -p $OUTDIR

echo "Starting main analysis workflow."

./gfa-to-fasta.py

./align-reads-to-edges.sh

./sam-to-bam.sh

# We need to have indexed the BAM file in order to use fetch() in Pysam
# in the partially-mapped reads filtering script below. (And, at least
# as far as I'm aware, we need to sort a BAM file before indexing it.)
./sort-and-index-bam.sh $OUTDIR/alignment.bam $OUTDIR/aln-sorted.bam

./filter-overlapping-supplementary-alignments.py

# I'm not sure if sorting this BAM again is strictly necessary, but we err
# on the side of safety and re-sort it anyway. We do need to index it to use
# pysam with it.
./sort-and-index-bam.sh $OUTDIR/overlap-supp-aln-filtered-aln.bam $OUTDIR/overlap-supp-aln-filtered-and-sorted-aln.bam

./filter-partially-mapped-reads.py

./sort-and-index-bam.sh $OUTDIR/pmread-filtered-aln.bam $OUTDIR/fully-filtered-and-sorted-aln.bam

./bam-to-pileup.py

echo "Done with the main analysis workflow!"
