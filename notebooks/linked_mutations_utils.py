import pileup

# This file defines various constants / functions used in the
# linked mutation analysis notebooks.

# unless (pos j) - (pos i) < this, we do not consider i and j linked.
MAX_DIST_BTWN_LINKED_POSITIONS_NONINCLUSIVE = float("inf")


# unless at least this many reads have mutations at both pos i and pos j,
# we do not consider i and j linked.
MIN_COV_OF_MUTATIONS_AT_LINKED_POSITIONS = 1000


# unless |Reads(i, -)| + |Reads(-, j)| < this fraction * |Reads(i, j)|,
# we do not consider i and j linked.
MAX_NONLINKED_MUTATED_FRACTION_NONINCLUSIVE = 0.2


# How we call a mutation: only if
# (# mismatches) / (# mismatches + # matches) > MINFREQ.
# Defaults to 0.5%.
MINFREQ = 0.005


# Used to initialize entries in the pospair2groupcts defaultdicts.
# This was originally a lambda function, but that breaks pickle:
# https://stackoverflow.com/a/16439720
# And it looks like it needs to be in scope when loading the pickled
# defaultdicts, anyway. So we use an ordinary function here instead.
def emptyListOf4():
    return [0, 0, 0, 0]


# Finds all mutated positions in a genome. Stored here to enable reuse.
def find_mutated_positions(seq):
    seq2pos2pileup = pileup.load()
    mutated_positions = []
    for pos, pcol in enumerate(seq2pos2pileup[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        # We can be strict and filter out positions that don't pass the
        # coverage filter for linked reads -- no sense including these.
        if cov >= MIN_COV_OF_MUTATIONS_AT_LINKED_POSITIONS:
            
            # Actually "call" mutations, the same way we do elsewhere in
            # these analyses (albeit maybe with different values of MINFREQ).
            # Of course, this isn't the only way to do this.
            if pileup.naively_call_mutation(pcol, MINFREQ):
                # We need 0-indexed positions so we can work with skbio /
                # pysam. We therefore decrease these by 1, since the pileup
                # positions are 1-indexed.
                mutated_positions.append(pos - 1)
    return mutated_positions
