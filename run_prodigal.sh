#! /usr/bin/env bash
# running with the default normal (in this version, "single") mode since these
# are long enough sequences. We're splitting them up into multiple runs per the
# Prodigal wiki's advice --
# https://github.com/hyattpd/prodigal/wiki/advice-by-input-type.
OUT=/Poppy/mfedarko/sheep_metagenome/prodigal_out/
prodigal -i 166.fasta  -o $OUT/166_genes.txt  -a $OUT/166_proteins.faa  -f sco
prodigal -i 6018.fasta -o $OUT/6018_genes.txt -a $OUT/6018_proteins.faa -f sco
prodigal -i 7998.fasta -o $OUT/7998_genes.txt -a $OUT/7998_proteins.faa -f sco
