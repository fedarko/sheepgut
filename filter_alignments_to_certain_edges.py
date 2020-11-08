#! /usr/bin/env python3
import re

IN_FILENAME = "/Poppy/mfedarko/sheep_metagenome/selected_alignment_less_soft_clipping.sam"
OUT_FILENAME = "/Poppy/mfedarko/sheep_metagenome/selected_alignment_less_sc_filtered.sam"

edges_we_care_about = ["edge_6018", "edge_7998", "edge_166"]

with open(IN_FILENAME, "r") as sam_file:
    out_text = ""
    i = 0
    for line in sam_file:
        if line.startswith("@"):
            out_text += line
            i += 1
        else:
            rname = line.split("\t")[2]
            for e in edges_we_care_about:
                if e in rname:
                    out_text += line
                    i += 1
                    break

        # We wanna avoid storing a lot of the multi-GB SAM file in memory, but
        # we also want to save time (it may be wasteful to write stuff out
        # every line). As a silly compromise, every 10,000 lines (that will be
        # included in the output), we write to the output file.
        if i > 10000:
            with open(OUT_FILENAME, "a") as out_file:
                out_file.write(out_text)
            out_text = ""
            i = 0

    # Take care of remaining stuff
    if i > 0:
        with open(OUT_FILENAME, "a") as out_file:
            out_file.write(out_text)
