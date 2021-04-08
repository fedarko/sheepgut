#! /usr/bin/env bash
# The purpose of running CheckM is to get a sense of edge sequence quality --
# e.g. are any really obvious chimeras, etc.

OUTDIR=/Poppy/mfedarko/sheepgut/data-processing-scripts/output

# Use of lineage_wf based on
# https://github.com/Ecogenomics/CheckM/wiki/Quick-Start
# -- may be overkill since I really just want QA, but ok with me.
#
# Also, use of -t 10 based on Misha's use of CheckM on the sheep gut dataset.
checkm lineage_wf -t 10 $OUTDIR/edge_fna/ $OUTDIR/checkm_out/
