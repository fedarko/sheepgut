#! /usr/bin/env bash
# running with the default normal (in this version, "single") mode since these
# are long enough sequences. We're splitting them up into multiple runs per the
# Prodigal wiki's advice --
# https://github.com/hyattpd/prodigal/wiki/advice-by-input-type.
OUT=/Poppy/mfedarko/sheepgut/data-processing-scripts/output
PROUT=$OUT/prodigal

prodigal -i $OUT/edge_1671.fasta -o $POUT/edge_1671_genes.sco -a $POUT/edge_1671_proteins.faa -f sco
prodigal -i $OUT/edge_2358.fasta -o $POUT/edge_2358_genes.sco -a $POUT/edge_2358_proteins.faa -f sco
prodigal -i $OUT/edge_6104.fasta -o $POUT/edge_6104_genes.sco -a $POUT/edge_6104_proteins.faa -f sco
