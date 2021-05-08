# This file defines various constants / functions used in the
# linked mutation analysis notebooks.

# unless (pos j) - (pos i) < this, we do not consider i and j linked.
MAX_DIST_BTWN_LINKED_POSITIONS_NONINCLUSIVE = 3000


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
# "seq2pos2matchct", etc are the same as they are in the mutation JSON
# notebook
def find_mutated_positions(seq, seq2pos2matchct, seq2pos2mismatchct):
    mutated_positions = []
    for pos in seq2pos2matchct[seq].keys():
        matchct = seq2pos2matchct[seq][pos]
        mismatchct = seq2pos2mismatchct[seq][pos]
        cov = mismatchct + matchct
        
        # We can be strict and filter out positions that don't pass the
        # coverage filter for linked reads -- no sense including these.
        if cov >= MIN_COV_OF_MUTATIONS_AT_LINKED_POSITIONS:
            
            # Actually "call" mutations, the same way we do elsewhere in
            # these analyses (albeit maybe with different values of MINFREQ).
            # Of course, this isn't the only way to do this.
            if (mismatchct / cov) > MINFREQ:
                mutated_positions.append(int(pos))
    return mutated_positions
