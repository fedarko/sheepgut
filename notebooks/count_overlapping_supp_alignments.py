import pysam
from itertools import combinations
from collections import defaultdict

bf = pysam.AlignmentFile("../main-workflow/output/fully-filtered-and-sorted-aln.bam", "rb")

SEQS = ["edge_6104", "edge_1671", "edge_2358"]
for seq in SEQS:
    print(f"Looking at {seq}...")
    read2refranges = defaultdict(list)
    read2atleast_one_supp_seen = defaultdict(bool)
    for ri, read in enumerate(bf.fetch(seq), 1):
        rn = read.query_name
        rng = range(read.reference_start, read.reference_end)
        read2refranges[rn].append(rng)
        if read.is_supplementary:
            read2atleast_one_supp_seen[rn] = True

    reads_with_supp_ct = 0
    overlaps = set()
    for r in read2refranges:
        if len(read2refranges[r]) > 1:
            reads_with_supp_ct += 1
            for combo in combinations(read2refranges[r], 2):
                if set(combo[0]) & set(combo[1]):
                    overlaps.add(r)

    numreads = len(read2refranges)
    print(f"{numreads:,} unique reads.")

    pctreadswithsupp = 100 * (reads_with_supp_ct / numreads)
    print(f"{reads_with_supp_ct:,} / {numreads:,} ({pctreadswithsupp:.2f}%) unique reads that have supp. alignments.")

    pctreadswithoverlaps = 100 * (len(overlaps) / numreads)
    print(f"{len(overlaps):,} / {numreads:,} ({pctreadswithoverlaps:.2f}%) unique reads that have supp. alignments where at least one pair of alignments overlaps.")

    for r in read2refranges:
        if len(read2refranges[r]) > 1:
            if not read2atleast_one_supp_seen[r]:
                print(f"Read {r} had no supps but multiple alignments???")
