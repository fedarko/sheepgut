#! /usr/bin/env python3
# Prove that there are 138 synonymous single-nucleotide mutations
# and 438 nonsynonymous single-nucleotide mutations from each of
# the unique 64 DNA codons. To justify a point we make at one point
# in the paper.

import skbio

codons = []
for i in "ACGT":
    for j in "ACGT":
        for k in "ACGT":
            codons.append(i + j + k)

si = 0
ni = 0
for c in codons:
    aa = str(skbio.DNA.translate(skbio.DNA(c)))
    n = 0
    for pos in (0, 1, 2):
        posnt = c[pos]
        for altnt in sorted(set("ACGT") - set(posnt)):
            # it should be possible to do this without checking pos and just
            # using a single fancy slice operation, but this is more foolproof
            # imo and i am nothing if not a fool a lot of the time
            if pos == 0:
                alt_codon = altnt + c[1:]
            elif pos == 1:
                alt_codon = c[0] + altnt + c[2]
            else:
                alt_codon = c[:2] + altnt

            aa2 = str(skbio.DNA.translate(skbio.DNA(alt_codon)))
            if aa2 == aa:
                si += 1
                print(f"{c} ({aa}) -> {alt_codon} ({aa2}) is syn")
            else:
                ni += 1
                print(f"{c} ({aa}) -> {alt_codon} ({aa2}) is nonsyn")
            n += 1
    if n != 9:
        # something went very wrong
        raise ValueError("each codon should only have 9 alt codons???")

print(f"sum of Si = {si}")
print(f"sum of Ni = {ni}")