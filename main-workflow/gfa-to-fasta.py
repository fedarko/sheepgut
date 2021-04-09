#! /usr/bin/env python3
# Converts a GFA assembly graph to a FASTA file of all sequences
# within the graph. Notably, this ignores connections between sequences
# in the graph (...those are considered elsewhere in the report -- we just
# care about the sequences).

# Identify input graph location -- it's configurable to make running this
# on another system less painful, hopefully
with open("../config/input-graph", "r") as graph_filename:
    IN_FILENAME = graph_filename.read().strip()

OUT_FILENAME = "output/all_edges.fasta"

print("Creating FASTA of all edges in the assembly graph...")

fasta_out = ""
with open(IN_FILENAME, "r") as gfa_file:
    for line in gfa_file:
        if line.startswith("S\t"):
            split = line.strip().split("\t")
            seq = split[2]
            fasta_out += ">{}\n".format(split[1])
            fasta_out += split[2] + "\n"

with open(OUT_FILENAME, "w") as fasta_file:
    fasta_file.write(fasta_out)

print("FASTA created.")
