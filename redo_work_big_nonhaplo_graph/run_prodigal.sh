#! /usr/bin/env bash
# running with the default normal (in this version, "single") mode since these
# are long enough sequences. We're splitting them up into multiple runs per the
# Prodigal wiki's advice --
# https://github.com/hyattpd/prodigal/wiki/advice-by-input-type.
OUT=/Poppy/mfedarko/sheep_metagenome/redo_work_big_nonhaplo_graph/prodigal_out
prodigal -i edge_1371.fasta -o $OUT/edge_1371_genes.sco -a $OUT/edge_1371_proteins.faa -f sco
prodigal -i edge_2358.fasta -o $OUT/edge_2358_genes.sco -a $OUT/edge_2358_proteins.faa -f sco
prodigal -i edge_6104.fasta -o $OUT/edge_6104_genes.sco -a $OUT/edge_6104_proteins.faa -f sco
