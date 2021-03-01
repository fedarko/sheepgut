#! /usr/bin/env python3
GRAPH_FILENAME = "/Poppy/mkolmogo/sheep_meta/flye_big_2.8/assembly_graph.gfa"

edges_to_find = [
    "10465", # 7998 in the old graph
    "2358", # 6018 in the old graph
    "1371", # 166 in the old graph
]

with open(IN_FILENAME, "r") as gfa_file:
    for line in gfa_file:
        if line.startswith("L\t"):
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
