# TODO use json loading notebook and gene utils notebook;
# adapt this into its own notebook
# that generates a misc-text file about number of ties in mutated positions
import json
import os

JSONPREFIX = "."

seq2pos2totalcov = {}
with open(os.path.join(JSONPREFIX, "seq2pos2totalcov.json"), "r") as jf:
    seq2pos2totalcov = json.load(jf)

seq2pos2matchct = {}
with open(os.path.join(JSONPREFIX, "seq2pos2matchct.json"), "r") as jf:
    seq2pos2matchct = json.load(jf)

seq2pos2mismatchct = {}
with open(os.path.join(JSONPREFIX, "seq2pos2mismatchct.json"), "r") as jf:
    seq2pos2mismatchct = json.load(jf)
    
seq2pos2mismatches = {}
with open(os.path.join(JSONPREFIX, "seq2pos2mismatches.json"), "r") as jf:
    seq2pos2mismatches = json.load(jf)
seq2pos2mismatches["edge_6104"]
g = [str(x) for x in range(1208927, 1210076)]
g
g[0]
g[-1]
seq2pos2mismatches["edge_6104"]
mm = [seq2pos2mismatches["edge_6104"][x] for x in g]
mm
for m in mm:
    for v in mm[m].values():
        if v > 100:
            print(m, mm[m])
for m in g:
    for v in seq2pos2mismatches["edge_6104"][m].values():
        if v > 100:
            print(m, seq2pos2mismatches["edge_6104"][m])
import numpy as np
np.argmax
help(_)
np.argmax({"G": 6, "T": 728})
help(np.argmax)
d = {"G": 6, "T": 728}
max(d, key=d.get)
help(dict.get)
def get_val(seq, pos, pseudo_variant_caller):
    """Does "variant calling" (note the quotation marks) given a seq, position, and variant caller function.
    
    See histogram_maker() for context on pseudo_variant_caller. This function was abstracted from that function
    to make doing this in different contexts easier.
    """ 
    mismatchct = seq2pos2mismatchct[seq][str(pos)]
    matchct = seq2pos2matchct[seq][str(pos)]
    # Note that, as mentioned above, this isn't the "true" coverage -- it can be zero
    # if e.g. all the reads covering a position are deletions
    cov = mismatchct + matchct
    if cov > 0:
        val = pseudo_variant_caller(cov, mismatchct)
    else:
        # Assign "val" of 0 - we do not have the data to detect if this position is a mutation.
        val = 0
    return val
pvc = lambda cov, mismatches: 1 if (mismatches / cov) > 0.005 else 0
0.005 * 100
for seq in SEQS:
    for pos in seq2pos2cov[seq]:
        if get_val(seq, pos, pvc):
            alts = seq2pos2mismatches[seq][pos]
            maxalt = max(alts, key=alts.get)
            maxalts = [a if alts[a] == alts[maxalt]]
            if len(maxalts) > 1:
                print(seq, pos, alts)
numties = 0
SEQS = ["edge_6104", "edge_1671", "edge_2358"]
for seq in SEQS:
    for pos in seq2pos2totalcov[seq]:
        if get_val(seq, pos, pvc):
            alts = seq2pos2mismatches[seq][pos]
            maxalt = max(alts, key=alts.get)
            maxalts = [a for a in alts if alts[a] == alts[maxalt]]
            if len(maxalts) > 1:
                print(seq, pos, alts)
                print(seq2pos2mismatchct[seq][pos])
                print(seq2pos2matchct[seq][pos])
                numties += 1
