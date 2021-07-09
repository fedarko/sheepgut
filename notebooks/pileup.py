# This file contains various utilities to help working with the
# seq2pos2pileup file created in ../main-workflow/bam-to-pileup.py.
#
# The motivation of this file is to store all of the pileup data in memory
# (...yeah, I know, I know) for a set of MAGs of interest. This prevents
# having to repeatedly mess around with pysam / look at FASTA files / etc.
#
# This could be compressed further (it's definitely not a good idea to use this
# with more than a few MAGs) -- but for just a handful of MAGs it's pretty
# useful.
#
# This file is structured as a dict mapping sequence names
# (e.g. "edge_6104") to a list of pileup entries for all 1-indexed positions
# within this sequence. The 0-th entry in each of these lists is None
# (I know, I know!) and subsequent entries describe the pileup stats at this
# position.
#
# Pileup entries are formatted like:
#
#   [
#     [A, C, G, T],
#     ri,
#     D
#   ]
#
# A, C, G, and T are integers indicating the number of aligned
# nucleotides (of A, C, G, and T, respectively) to this reference position.
#
# ri is an integer in the range [0, 3]. This indicates the reference nucleotide
# at this position in the sequence. 0 -> A, 1 -> C, 2 -> G, 3 -> T. (This can
# be used to index into the [A, C, G, T] list of the pileup entry.)
#
# D indicates the number of aligned deletions to this position (in case we need
# it).
#
# Although it's a very simple format, this makes it easy to repeatedly do
# naive variant calling / plot mutation spectra / etc.


import pickle


def load(picklepath="../main-workflow/output/seq2pos2pileup.pickle"):
    # https://stackoverflow.com/a/18261955
    with open(picklepath, "rb") as picklefile:
        return pickle.load(picklefile)


def get_mismatch_cts(pileup):
    """Returns pileup[0] with the pileup[1]-th element removed.

    e.g. if pileup[0] = [99, 0, 30, 14] and pileup[1] = 2,
    this'll return [99, 0, 14].

    This corresponds to filtering the list of [A, C, G, T] aligned to a
    position to remove whichever of the four nucleotides corresponds to
    the reference nucleotide at this position.
    """
    ref_idx = pileup[1]
    return pileup[0][:ref_idx] + pileup[0][ref_idx + 1:]


def get_cov(pileup, raise_error_if_0x=False):
    """Returns the sum of matches + mismatches for a position.

    Optionally raises an error if this particular pileup has a 
    (match + mismatch) coverage of 0.
    """
    cov = sum(pileup[0])
    if raise_error_if_0x and cov <= 0:
        raise ValueError(f"pileup {pileup} has coverage of {cov}x.")
    else:
        return cov


def get_max_freq_alt_nt(pileup):
    """Raises an error if there are no mismatches.

    Breaks ties arbitrarily.
    """
    ref_idx = pileup[1]

    alts = {}
    max_alt_freq = 0
    for ni, nt in enumerate("ACGT"):
        if ni != ref_idx:
            alts[nt] = pileup[0][ni]
            max_alt_freq = max(max_alt_freq, alts[nt])   

    if max_alt_freq == 0:
        raise ValueError("No mismatches at this position in the pileup.")

    # Retrieve max-freq alternate nucleotide. Based on
    # https://stackoverflow.com/a/280156.
    # (Note that if there's a tie, the result is arbitrary. Shouldn't 
    # be a big deal.)
    max_freq_alt_nt = max(alts, key=alts.get)  

    # Warn if we need to arbitrarily break a tie
    if list(alts.values()).count(max_alt_freq) > 1:
        print(f"Multiple max-freq alt nucleotides: {alts} for pileup {pileup}")
        print(f"\t(Arbitrarily breaking tie: selecting max alt = {max_freq_alt_nt}.)")

    return max_freq_alt_nt
    

def get_mismatch_pcts(pileup):
    """Like get_mismatch_cts(), but returns percentages.

    The percentages include all aligned nucleotides (both matches and
    mismatches).
    """
    cov = get_cov(pileup, raise_error_if_0x=True)
    return [c / cov for c in get_mismatch_cts(pileup)]


def get_agg_mismatch_pct(pileup):
    """Returns the total percentage of mismatched reads aligned to a pos.

    So, e.g., if the reference nt is "A" and the pileup counts are
    [100, 20, 30, 50], then this'll return (20+30+50) / (100+20+30+50)
    = (100 / 200) = 0.5.

    THIS PROBABLY SHOULDN'T BE USED FOR VARIANT CALLING YO because it
    aggregates all the alternate nucleotides at a position. I mean, you could
    totally call variants that way, but it disagrees with what we describe in
    our paper (see: definition of "freq(pos)").
    """
    return sum(get_mismatch_pcts(pileup))

def get_max_mismatch_pct(pileup):
    """Returns the *maximum* percentage of an alternate nt aligned to a pos.

    So, e.g., if the reference nt is "A" and the pileup counts are
    [100, 20, 30, 50], then this'll return (50) / (100+20+30+50)
    = (50 / 200) = 0.25.

    Counterpart to get_agg_mismatch_pct(). This is equivalent to freq(pos),
    as of writing.
    """
    return max(get_mismatch_pcts(pileup))


def get_max_freq_alt_nt_pct(pileup):
    # Fancier version of naively_call_mutation().
    # Useful for checking lots of different values of p
    # for the same position.
    mismatch_cts = get_mismatch_cts(pileup)
    max_freq_alt_nt_ct = max(mismatch_cts)
    cov = get_cov(pileup)
    if cov > 0:
        return (max_freq_alt_nt_ct / cov)
    else:
        # Arguably this should be undefined, since it's n / 0 -- but
        # the desired result for these situations (we don't call a mutation
        # because this position is completely uncovered) is respected by
        # just treating the max-freq alternate nucleotide percentage as 0.
        return 0
    

def naively_call_mutation(pileup, p):
    # p should be in [0, 1].
    max_freq_alt_nt_pct = get_max_freq_alt_nt_pct(pileup)
    return max_freq_alt_nt_pct > p
