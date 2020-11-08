#! /usr/bin/env bash

# TO RUN THIS SCRIPT: Modify the edges_to_keep variable in
# convert_gfa_to_fasta.py to select just the edge IDs you want to keep. 
# Then just run this using ./produce_spectra.sh.

./convert_gfa_to_fasta.py
echo "Did GFA->FASTA conversion."
./align_reads_to_edges.sh
echo "Did alignment."
./remove_softclipping_in_sam.py
echo "Filtered to fully kept reads."
./view_sort_index.sh
echo "Did v/s/i."
./compute_mutations.sh
echo "Done."
