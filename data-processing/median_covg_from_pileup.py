#! /usr/bin/env python3

from collections import defaultdict
from statistics import median

seq2poscovs = defaultdict(list)
with open("pileup.txt", "r") as pf:
    for line in pf:
        s = line.strip().split()
        seq2poscovs[s[0]].append(int(s[3]))

for seq in seq2poscovs.keys():
    print("Seq {}: Median coverage of {}.".format(
        seq, median(seq2poscovs[seq])
    ))
