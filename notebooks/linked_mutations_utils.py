import pileup
from collections import defaultdict

# This file defines various constants / functions used in the
# linked mutation analysis notebooks.

# We only consider a position if it has at least this much coverage. This is a
# very greedy thing to do, but thanks to super deep coverage in these genomes
# it's feasible.
MINCOV = 1000

# We only add an allele node to the graph if its frequency is greater than
# this. Note that this is >, not >= -- so a value of 1 here means that this
# must be at least 2 (i.e. this nt at this position in this sequence was
# observed at least twice).
MIN_ALLELE_FREQ_EXCLUSIVE = 1

# In order for us to connect two allele nodes, at least this many reads must
# span both allele's positions (this isn't the only condition for creating an
# edge; see MINLINK_EXCLUSIVE below).
MINSPAN = 500

# In order for us to connect two allele nodes, link(i, j, N_i, N_j) between
# positions i and j with nucleotides N_i and N_j must be GREATER THAN this (we
# use > instead of >= since this can be 0).
MINLINK_EXCLUSIVE = 0

# Value of p used for calling a position as "mutated" or not. defaults to 0.5%
p = 0.5 / 100

def gen_ddi():
    """Returns a new defaultdict(int).

    Needed because pickle can't handle lambda functions.
    """
    return defaultdict(int)

# Finds all mutated positions in a genome. Stored here to enable reuse.
def find_mutated_positions(seq):
    """Returns a list of mutated positions in a genome.

    Mutated positions are stored in the list as 0-indexed integers (so the
    first position in a genome, if included, would be 0).

    Mutated positions are classified as being mutated based on
    pileup.naively_call_mutation() with the value of p listed above. We also
    only include mutated positions that have coverage of at least MINCOV.
    """
    seq2pos2pileup = pileup.load()
    mutated_positions = []
    for pos, pcol in enumerate(seq2pos2pileup[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        # We can be strict and filter out positions that don't pass the
        # coverage filter for linked reads -- no sense including these.
        if cov >= MINCOV:
            if pileup.naively_call_mutation(pcol, p):
                # We need 0-indexed positions so we can work with skbio /
                # pysam. We therefore decrease these by 1, since the pileup
                # positions are 1-indexed.
                mutated_positions.append(pos - 1)
    return mutated_positions
