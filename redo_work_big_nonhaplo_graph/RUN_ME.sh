#! /usr/bin/env bash

#./convert_gfa_to_fasta.py
#echo "Created FASTA of all edges in the assembly graph."
echo "Aligning reads to edges..."
./align_reads_to_edges.sh
echo "Did alignment."
./filter_alignment.sh
echo "Did filtering."
#./remove_softclipping_in_sam.py
#echo "Filtered to fully (ish) aligned reads."
./view_sort_index.sh
echo "Did v/s/i."
./bam2pointinfojsons.py
#./compute_mutations.sh
#echo "Computed pileup."
#./pileup2pointinfo.py
./zip_jsons.sh
echo "Done."
