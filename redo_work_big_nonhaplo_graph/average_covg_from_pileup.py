#! /usr/bin/env python3

from collections import defaultdict

# Compute "real" coverages -- the average number of reads aligned to each
# position in each sequence.
# Defaultdicts of int have 0 as their default value, which is what we want
seq2totalnumalignedreadsateachposition = defaultdict(int)
seq2len = defaultdict(int)
with open("pileup.txt", "r") as pf:
    for line in pf:
        s = line.strip().split()
        seq2totalnumalignedreadsateachposition[s[0]] += int(s[3])
        seq2len[s[0]] += 1

for seq in seq2len.keys():
    print("Seq {}: Avg coverage of {}.".format(
        seq, seq2totalnumalignedreadsateachposition[seq] / seq2len[seq]
    ))
