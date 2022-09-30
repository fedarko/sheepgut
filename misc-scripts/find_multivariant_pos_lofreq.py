#! /usr/bin/env python3
import pysam
from collections import Counter

# maps seq name --> list of all mutations' positions. includes duplicates
# if multiple mutations (e.g. --> A, --> C) occur at a given position.
seq2seenpos = {"edge_6104": [], "edge_1671": [], "edge_2358": []}

# seq name --> pos --> list of alt nts
seq2pos2alts = {"edge_6104": {}, "edge_1671": {}, "edge_2358": {}}
# seq name --> pos --> ref
seq2pos2ref = {"edge_6104": {}, "edge_1671": {}, "edge_2358": {}}

vf = pysam.VariantFile("../seqs/lofreq.vcf")

# I'm pretty sure the positions in a vcf file are usually in sorted order, so
# we could use that to make this code much more efficient. however, it already
# takes like a few seconds max to run, so ...not gonna bother with that r/n
for r in vf:
    seq2seenpos[r.contig].append(r.pos)
    alt = r.alts[0].upper()
    if r.pos in seq2pos2alts[r.contig]:
        if alt in seq2pos2alts[r.contig][r.pos]:
            raise ValueError(
                "Weird: completely duplicate mutation (w/ same alt nt) at pos "
                f"{r.pos} in contig {r.contig}"
            )
        else:
            seq2pos2alts[r.contig][r.pos].append(alt)
            if r.ref.upper() != seq2pos2ref[r.contig][r.pos]:
                raise ValueError(
                    f"Weird: inconsistent ref nt at pos {r.pos} in contig "
                    f"{r.contig}"
                )
    else:
        seq2pos2alts[r.contig][r.pos] = [alt]
        seq2pos2ref[r.contig][r.pos] = r.ref.upper()

for seq in seq2seenpos:
    hl = "=" * 79
    print(f"{hl}\nSequence {seq}\n{hl}")
    num_multimuts = 0
    mut_pos_ctr = Counter(seq2seenpos[seq])
    for mut_pos in mut_pos_ctr.keys():
        mc = mut_pos_ctr[mut_pos]
        if mc > 1:
            ref = seq2pos2ref[seq][mut_pos]
            alts = seq2pos2alts[seq][mut_pos]
            print(f"\tPosition {mut_pos:,} has {mc} alt nts: {ref} -> {alts}")
            num_multimuts += 1
    print("\t" + ("-" * 45))
    print(f"\t{seq} has {num_multimuts:,} multi-mutation position(s).")
