#!/usr/bin/env python3

from collections import defaultdict

FASTALOC = "/Poppy/mkolmogo/sheep_meta/flye_big_2.8/assembly.fasta"
OUTLOC = "/Poppy/mfedarko/sheep_metagenome/redo_work_big_nonhaplo_graph/newseqs.fasta"

output_str = ""

names = ["scaffold_10465", "contig_2358", "contig_1371"]

in_selected_seqs = False
with open(FASTALOC, "r") as ff:
    for line in ff:
        if line.startswith(">"):
            declaration_is_match = False
            for n in names:
                if line.strip() == ">{}".format(n):
                    declaration_is_match = True
                    in_selected_seqs = True
                    break
            if not declaration_is_match:
                in_selected_seqs = False
        if in_selected_seqs:
            output_str += line

with open(OUTLOC, "w") as of:
    of.write(output_str)
