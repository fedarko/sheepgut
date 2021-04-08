#! /usr/bin/env python3
# TODO: don't actually do this manually. just align graph edges' fasta to
# 166/6018/7998 fasta. oh can even use assembly.fasta LOLLLLL
GRAPH_FILENAME = "/Poppy/mkolmogo/sheep_meta/flye_big_2.8_haplo/assembly_graph.gfa"
FASTA_FILENAME = "/Poppy/mfedarko/sheep_metagenome/166_6018_7998.fasta"

edges_to_find = [
    "7998",
    "166",
    "6018"
]

for e in edges_to_find:

with open(IN_FILENAME, "r") as gfa_file:
    for line in gfa_file:
        if line.startswith("S\t"):
            for e in edges_to_find:
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
