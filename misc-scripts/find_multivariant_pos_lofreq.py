#! /usr/bin/env python3
import pysam
from collections import defaultdict, Counter

# maps seq name --> list of all mutations' positions. includes duplicates
# if multiple mutations (e.g. --> A, --> C) occur at a given position.
seq2seenpos = defaultdict(list)

vf = pysam.VariantFile("../seqs/lofreq.vcf")

# I'm pretty sure the positions in a vcf file are usually in sorted order, so
# we could use that to make this code much more efficient. however, it already
# takes like a few seconds max to run, so ...not gonna bother with that r/n
for r in vf:
    seq2seenpos[r.contig].append(r.pos)

for seq in seq2seenpos:
    hl = "=" * 79
    print(f"{hl}\nSequence {seq}\n{hl}")
    num_multimuts = 0
    mut_pos_ctr = Counter(seq2seenpos[seq])
    for mut_pos in mut_pos_ctr.keys():
        mc = mut_pos_ctr[mut_pos]
        if mc > 1:
            print(f"\tPosition {mut_pos:,} has {mc} mutations")
            num_multimuts += 1
    print("\t" + ("-" * 45))
    print(f"\t{seq} has {num_multimuts:,} multi-mutation positions.")
