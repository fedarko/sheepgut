#! /usr/bin/env bash
# The purpose of running CheckM is to get a sense of edge sequence quality --
# e.g. are any really obvious chimeras, etc.

OUTDIR=/Poppy/mfedarko/sheepgut/data-processing-scripts/output

# Use of lineage_wf based on
# https://github.com/Ecogenomics/CheckM/wiki/Quick-Start
# -- may be overkill since I really just want QA, but ok with me.
#
# Also, use of -t 10 based on Misha's use of CheckM on the sheep gut dataset.
#
# edge_fna/ should be a directory containing the sequences of interest; these
# seqs should have .fna extensions. (The reason we don't just tell CheckM to
# look at all FASTA files in the output directory is that there could be lots
# of other stuff in there...)
checkm lineage_wf -t 10 $OUTDIR/edge_fna/ $OUTDIR/checkm/
