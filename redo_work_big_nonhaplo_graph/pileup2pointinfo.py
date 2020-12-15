#! /usr/bin/env python3
# Converts a pileup file produced by "samtools mpileup" to JSON files.
# These JSON files are much smaller (so easier to download from the server) and
# reduce the need to continuously run through the pileup file every time I want
# to rerun an analysis using this data.

import re
import json
from collections import defaultdict
from parse_pileup_read_bases import parse_read_bases

# Hard-coded stuff. Ideally this'd be extracted from a FASTA file or something.
SEQ2LEN = {
    "edge_1371": 1634973,
    "edge_2358": 2806161,
    "edge_6104": 1289244,
}

seq2pos2cov = defaultdict(dict)
seq2pos2matches = defaultdict(dict)
seq2pos2nonmatches = defaultdict(dict)

with open("pileup.txt", "r") as pf:
    i = 0
    prev_pos = None
    prev_seq = None
    for line in pf:
        s = line.strip().split("\t")
        seq = s[0]
        pos = int(s[1])
        # Try to keep track of discontinuities in the pileup. If stuff gets
        # skipped, be liberal with errors.
        if prev_pos is None:
            prev_pos = pos
            prev_seq = seq
        else:
            if pos != prev_pos + 1:
                if seq == prev_seq:
                    raise ValueError(
                        "Position just skipped from {} to {} within seq {}.".format(
                            prev_pos, pos, seq
                        )
                    )
                else:
                    if pos == 1:
                        if prev_pos == SEQ2LEN[prev_seq]:
                            # We just finished going over one sequence and now
                            # we're starting a new one. That's fine.
                            prev_pos = pos
                            prev_seq = seq
                        else:
                            raise ValueError("Exited seq {} early.".format(prev_seq))
                    else:
                        raise ValueError(
                            "Skipped from pos {} in {} to pos {} in {}?".format(
                                prev_pos, prev_seq, pos, seq
                            )
                        )
            else:
                # We're at the previous position + 1, so we are probably fine.
                # However, a sneaky extra case is that we just jumped to the
                # expected position but at a completely different sequence.
                # That's no good!
                if seq != prev_seq:
                    raise ValueError(
                        (
                            "Position incremented as expected from {} to {}, but "
                            "seq changed from {} to {}???"
                        ).format(prev_pos, pos, prev_seq, seq)
                    )
                # If we made it *here*, then pos == prev_pos + 1 and
                # seq == prev_seq. So, we're good.

        # ... Now that things seem good, let's record this line.
        cov = int(s[3])
        seq2pos2cov[seq][pos] = cov
        if cov > 0:
            matchct, non_matches = parse_read_bases(s[4])
            seq2pos2matches[seq][pos] = matchct
            seq2pos2nonmatches[seq][pos] = non_matches
        else:
            # Right now, I expect that all sequences should be covered
            # throughout. If this isn't the case, report it.
            # (Note that uncovered positions may also just straight up not be
            # listed in the pileup file, which is something I've seen before.
            # As a TODO we could detect this by checking to see if a position
            # "skips.")
            raise ValueError(
                "Found uncovered position: pos {} in {}".format(pos, seq)
            )
        i += 1
        if i % 500000 == 0:
            print("On line {}".format(i))

# yes i know this is extremely lazy
with open("seq2pos2cov.json", "w") as jf:
    jf.write(json.dumps(seq2pos2cov))

with open("seq2pos2matches.json", "w") as jf:
    jf.write(json.dumps(seq2pos2matches))

with open("seq2pos2nonmatches.json", "w") as jf:
    jf.write(json.dumps(seq2pos2nonmatches))
