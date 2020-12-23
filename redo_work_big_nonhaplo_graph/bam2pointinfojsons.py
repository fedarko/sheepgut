#! /usr/bin/env python3
# Converts a BAM file to JSON files that describe the pileup.
# These JSON files are much smaller than either the plaintext pileup output
# produced by "samtools mpileup" or the original BAM file, so they're a lot
# easier to download from the server (as well as to read).

import json
from collections import defaultdict
import pysam
import pysamstats

# Hard-coded stuff. Ideally this'd be extracted from a FASTA file or something.
SEQ2LEN = {
    "edge_1371": 1634973,
    "edge_2358": 2806161,
    "edge_6104": 1289244,
}

# Coverage from the alignment -- includes ALL stuff from the alignment, except
# for stuff we manually filtered before producing the BAM file (i.e.
# partially-mapped reads).
#
# So this means that insertions / deletions are included here. We COULD just
# use matches + mismatches to determine this, but I don't think that makes
# sense...?
# This does mean that we can no longer infer mismatches from just (coverage)
# and (match count), hence why we have added an additional dict for
# seq2pos2mismatchct. I think this is a reasonable way of handling this.
seq2pos2totalcov = defaultdict(dict)

# The number of matches aligned to a position.
seq2pos2matchct = defaultdict(dict)

# The number of mismatches aligned to a position.
# Should be equal to the sum of seq2pos2nonmatches.values() at this seq & pos
# -- this is just done here for convenience.
seq2pos2mismatchct = defaultdict(dict)

# Another dict, mapping mismatched nucleotides to their frequency in the
# alignment at a position. 0s are omitted for the sake of filesize.
seq2pos2mismatches = defaultdict(dict)

bf = pysam.AlignmentFile("aln-sorted.bam")

for seq in SEQ2LEN.keys():
    # We use start=0 and end=(sequence length) because the start and end params
    # of pysamstats.stat_variation() are 0-indexed (although the normal
    # samtools pileup format's coordinates are 1-indexed).
    # Also, we set max_depth at 100k because having it low enough silently
    # limits coverage, I guess???? asodfij. See pysamstats docs.
    for i, rec in enumerate(pysamstats.stat_variation(
        bf, chrom=seq, fafile="newseqs.fasta", start=0, end=SEQ2LEN[seq],
        truncate=True, max_depth=100000
    ), 1):
        if rec["N"] > 0:
            raise ValueError("Hang on, there shouldn't be any Ns in this data")

        # Output stuff in the JSONs as 1-indexed positions (because that's what
        # my analysis code assumes, because that's what mpileup produces)
        if rec["pos"] + 1 != i:
            raise ValueError(
                "Found discontinuity in traversal: {}-th pos, but rec['pos'] + 1 is {}".format(
                    i, rec["pos"] + 1
                )
            )
        pos = i

        seq2pos2totalcov[seq][pos] = rec["reads_all"]
        seq2pos2matchct[seq][pos] = rec["matches"]
        seq2pos2mismatchct[seq][pos] = rec["mismatches"]

        non_matches = {}
        possible_non_matches = set("ACGT") - set(rec["ref"])
        for poss_non_match in possible_non_matches:
            if rec[poss_non_match] > 0:
                non_matches[poss_non_match] = rec[poss_non_match]
        seq2pos2mismatches[seq][pos] = non_matches

        if i % 100000 == 0:
            print("Just processed {} positions in {}".format(i, seq))

    print("Just processed seq {}".format(seq))

# ...yes i know this is extremely lazy
with open("seq2pos2totalcov.json", "w") as jf:
    jf.write(json.dumps(seq2pos2totalcov))

with open("seq2pos2matchct.json", "w") as jf:
    jf.write(json.dumps(seq2pos2matchct))

with open("seq2pos2mismatchct.json", "w") as jf:
    jf.write(json.dumps(seq2pos2mismatchct))

with open("seq2pos2mismatches.json", "w") as jf:
    jf.write(json.dumps(seq2pos2mismatches))
