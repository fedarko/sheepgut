#! /usr/bin/env python3
IN_FILENAME = "/Poppy/mkolmogo/sheep_meta/flye_2.8_haplo/assembly_graph.gfa"
OUT_FILENAME = "/Poppy/mfedarko/sheep_metagenome/selected_edges.fasta"

# TODO: either take in these from a file, or load using argparse
edges_to_keep = [
    "7998",
    "166",
    "6018", "14948", "1809", "6017", "6021", "6019",
    "2069",
    "1829", "1879", "3735", "3736", "8360", "17289", "22832", "22833", "22834", "23160", "24017", "25923", "26447", "35242", "35953"
]

fasta_out = ""
with open(IN_FILENAME, "r") as gfa_file:
    for line in gfa_file:
        if line.startswith("S\t"):
            is_selected = False
            for e in edges_to_keep:
                if "edge_{}\t".format(e) in line:
                    is_selected = True
                    break
            if is_selected:
                split = line.strip().split("\t")
                seq = split[2]
                fasta_out += ">{}_len{}_cov{}\n".format(
                    split[1], len(seq), split[3][5:]
                )
                fasta_out += split[2] + "\n"

with open(OUT_FILENAME, "w") as fasta_file:
    fasta_file.write(fasta_out)
