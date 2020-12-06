#! /usr/bin/env python3
# Converts a pileup file produced by "samtools mpileup" to JSON files.
# The goal of this is
# making the analysis of point mutations easier (i.e. we don't care about
# indels).
#
# UPDATE: actually just dumping to json was easier lmao, ignore the rest of
# these docs
#
# Each TSV file produced has three columns:
# position (tab) num of matches aligned (tab) total num of reads aligned
#
# position: This will be 1-indexed, since the pileup format is 1-indexed (see
# http://www.htslib.org/doc/samtools-mpileup.html#Pileup_Format).
#
# num of matches aligned: This is equal to the number of "." and "," characters
# present within the "Read bases" column in the pileup file at this position.
# Again, indels do not seem to impact the number of "." and ","s, so we can
# just count these literally.
#
# total num of reads aligned: This is just equal to the "Number of reads
# covering this position" column in the pileup file at this position. This
# number seems to not include indels, so that's nice -- that's what we want.

import json
from collections import defaultdict

seq2pos2cov = defaultdict(dict)
seq2pos2matches = defaultdict(dict)

with open("pileup.txt", "r") as pf:
    i = 0
    for line in pf:
        s = line.strip().split("\t")
        num_matches = s[4].count(",") + s[4].count(".")
        pos = int(s[1])
        seq = s[0]
        seq2pos2cov[seq][pos] = int(s[3])
        seq2pos2matches[seq][pos] = num_matches
        i += 1
        if i % 500000 == 0:
            print("On line {}".format(i))

with open("seq2pos2cov.tsv", "w") as jf:
    jf.write(json.dumps(seq2pos2cov))

with open("seq2pos2matches.tsv", "w") as jf:
    jf.write(json.dumps(seq2pos2matches))
