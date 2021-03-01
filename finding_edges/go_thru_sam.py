#! /usr/bin/env python3

import re
from collections import defaultdict

matchhit = re.compile("(\d+)M")

def total_matches_in_cigar(cigar):
    matches = matchhit.findall(cigar)
    if matches:
        total_matches = sum([int(c) for c in matches])
        return total_matches
    return 0

seq2matches = defaultdict(list)
with open("graph_edges_to_166_etal_nonhaplo.sam", "r") as samfile:
    i = 0
    for line in samfile:
        if not line.startswith("@"):
            s = line.split("\t")
            for seq in ("166", "6018", "7998"):
                if "edge_{}".format(seq) in s[2]:
                    if total_matches_in_cigar(s[5]) > 100000:
                        seq2matches[seq].append((s[0], s[5]))
        i += 1
        if i % 1000 == 0:
            print("processed {} lines".format(i))

for seq in ("166", "6018", "7998"):
    print(seq, seq2matches[seq])
