#! /usr/bin/env python3
import re

IN_FILENAME = "/Poppy/mfedarko/sheep_metagenome/5xalignment.sam"
OUT_FILENAME = "/Poppy/mfedarko/sheep_metagenome/5xalignment_less_soft_clipping.sam"

softclipping = re.compile("(\d+)S")
with open(IN_FILENAME, "r") as sam_file:
    out_text = ""
    i = 0
    for line in sam_file:
        # ignore the header lines
        if not line.startswith("@"):
            cigar = line.split("\t")[5]
            matches = softclipping.findall(cigar)
            if matches:
                total_softclipped_posns = sum([int(c) for c in matches])
                if total_softclipped_posns < 50:
                    out_text += line
                    i += 1
            else:
                out_text += line
                i += 1
        else:
            out_text += line
            i += 1

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
