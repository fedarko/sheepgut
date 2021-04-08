#! /usr/bin/env python3
# Simple code to find weird stuff (e.g. positions that seem to be uncovered) in
# the pileup.
# Horribly unoptimized, but convenient for quick validation.

#with open("pileup.txt", "r") as pf:
#    for line in pf:
#        if line.startswith("edge_6104"):
#            if "208" in line:
#                with open("e6104p208.txt", "w") as ff:
#                    ff.write(line)
#                    break

with open("pileup.txt", "r") as pf:
    stuff_to_save = ""
    i = 0
    saving = False
    for line in pf:
        if line.startswith("edge_6104"):
            if "\t200" in line:
                saving = True
            if saving:
                stuff_to_save += line
                i += 1
                if i > 200:
                    break
    with open("saved-stuff.txt", "w") as sf:
        sf.write(stuff_to_save)
