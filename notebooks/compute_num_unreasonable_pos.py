import pileup
from collections import defaultdict


SEQS = ["edge_6104", "edge_1671", "edge_2358"]
seq2pos2pileup = pileup.load()
seq2unreasonable = defaultdict(int)
for seq in SEQS:
    print(f"On {seq}...")
    for pos, pospileup in enumerate(seq2pos2pileup[seq][1:], 1):
        if not pileup.is_reasonable(pospileup):
            seq2unreasonable[seq] += 1
            # print(f"{seq}: pos {pos:,} unreasonable")

print("=" * 79)
print("Totals")
print("=" * 79)
for seq in SEQS:
    print(f"{seq}: {seq2unreasonable[seq]:,} unreasonable positions")
