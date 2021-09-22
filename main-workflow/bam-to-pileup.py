#! /usr/bin/env python3
# Converts a BAM file to a simple, gross file describing the pileup.
# See notebooks/pileup.py in this repo for a description of how this file
# is structured.

import pickle
from collections import defaultdict
import skbio
import pysam
import pysamstats

# Script configuration settings:
# 1. Which seqs to use? these should have FASTA files in ../seqs/ (could alter
#    this script to instead look at the all_edges.fasta file in output/, but
#    that might be slow if done a bunch of times)
# 2. How often (i.e. after processing how many positions of a given sequence)
#    should we update the user about what's happening?
SEQS = ["edge_6104", "edge_1671", "edge_2358"]
UPDATE_FREQ = 5000

print("Converting the BAM file to a pickled dict...")

print(f"Using seqs {SEQS}.")

seq2pos2pileup = defaultdict(dict)

bf = pysam.AlignmentFile("output/fully-filtered-and-sorted-aln.bam")

# Maps ACGT to their position in ["A", "C", "G", "T"] -- used to easily
# separate the # of reference matches from the # of other matches
nt2idx = {"A": 0, "C": 1, "G": 2, "T": 3}

for seq in SEQS:
    fasta = str(skbio.DNA.read(f"../seqs/{seq}.fasta"))
    seqlen = len(fasta)
    print(f"Seq {seq} has length {seqlen:,} bp.")

    # We use 1-indexing in our output pileup data, even though python lists are
    # 0-indexed -- so the leading None lets us do this easily (I know that
    # doing this is horrible, but there's a lot of code I have which assumes
    # 1-indexing... Will see if I have time to fix this at some point in the
    # future. Maybe in a few decades.)
    pos2pileup = [None]

    # We use start=0 and end=(sequence length) because the start and end params
    # of pysamstats.stat_variation() are 0-indexed (although the normal
    # samtools pileup format's coordinates are 1-indexed, and although our
    # output will be 1-indexed).
    # Also, we set max_depth at 1 mil because having it low enough silently
    # limits coverage, I guess???? asodfij. See pysamstats docs.
    # Lastly, we use pad=True so that even uncovered positions are included --
    # edge 1671 has some regions in the middle that are uncovered after read
    # filtering, at least as of writing.
    for pos, rec in enumerate(pysamstats.stat_variation(
        bf, chrom=seq, fafile="output/all_edges.fasta", start=0, end=seqlen,
        truncate=True, max_depth=1000000, pad=True
    ), 1):
        if rec["N"] > 0:
            raise ValueError("Hang on, there shouldn't be any Ns in this data")

        rpos = rec["pos"] + 1
        if rpos != pos:
            raise ValueError(
                f"Found discontinuity in traversal: {pos}-th pos, but "
                f"rec['pos'] + 1 is {rpos}"
            )

        matches = rec["matches"]
        mismatches = rec["mismatches"]

        ref_nt = rec["ref"]
        if ref_nt != fasta[pos - 1]:
            raise ValueError(
                "Looks like FASTA and pysamstats disagree on ref nt at pos "
                f"{pos - 1} (0-indexed), a.k.a. {pos} (1-indexed). "
                f"FASTA says it's {fasta[pos-1]}; pysamstats says it's "
                f"{ref_nt}."
            )

        ref_idx = nt2idx[ref_nt]

        pos2pileup.append([
            [rec["A"], rec["C"], rec["G"], rec["T"]],
            ref_idx,
            rec["deletions"]
        ])

        # Print occasional status updates for my own sanity
        if pos % UPDATE_FREQ == 0:
            pct = 100 * (pos / seqlen)
            print(
                f"Just processed {pos:,} / {seqlen:,} ({pct:.2f}%) "
                f"positions in {seq}."
            )
            print(f"Pileup of position {pos:,} is:\n\t{pos2pileup[pos]}")

    seq2pos2pileup[seq] = pos2pileup
    print(f"Just processed seq {seq}.")

with open("output/seq2pos2pileup.pickle", "wb") as pf:
    pf.write(pickle.dumps(seq2pos2pileup))

print("Created pickle file.")
