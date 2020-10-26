#! /usr/bin/env bash

# TO RUN THIS SCRIPT: Modify the edges_to_keep variable in
# convert_gfa_to_fasta.py to select just the edge IDs you want to keep. 
# Then just run this using ./produce_spectra.sh.

./convert_gfa_to_fasta.py
./align_reads_to_edges.sh
./view_sort_index.sh
./compute_mutations.sh
