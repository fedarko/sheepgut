#! /usr/bin/env python3
IN_FILENAME = "/Poppy/mkolmogo/sheep_meta/flye_big_2.8/assembly_graph.gfa"
OUT_FILENAME = "output/all_edges.fasta"

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
