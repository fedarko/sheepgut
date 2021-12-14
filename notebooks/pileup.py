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
from pleuk_copied_code import get_alt_pos_info


# We consider a p-mutation with freq(pos) >= this percentage to be
# "high-frequency", aka indisputable. Depending on what we are using mutation
# calling for, we may want to not call high-frequency mutated positions (for
# example, this will throw off the computation of decoy genomes).
#
# This should be a value in the range (0, 50], analogous to the "p" value that
# is a parameter for naively_call_mutation() below.
HIGH_FREQUENCY_MIN_PCT = 5


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
    (match + mismatch) coverage of 0. (The error will also come up
    if through some twist of fate this pileup has a coverage of less
    than 0. If that happens, ...try calling an exorcist?)
    """
    cov = sum(pileup[0])
    if raise_error_if_0x and cov <= 0:
        raise ValueError(f"pileup {pileup} has coverage of {cov}x.")
    else:
        return cov


def get_alt_info_from_pleuk(pileup, warn_if_tie=False):
    """Convenience function -- takes a pileup and calls Pleuk's
    get_alt_pos_info() function on it. See the Pleuk bam_utils library
    for more details about this."""

    cts = pileup[0]
    nt2ct = {"A": cts[0], "C": cts[1], "G": cts[2], "T": cts[3]}
    cov, alt_freq, alt_nt = get_alt_pos_info(nt2ct)

    if warn_if_tie:
        if cts.count(alt_freq) > 1:
            print(f"Warning about pileup {pileup}:")
            print(
                "\tMultiple nucleotides w/ same freq as second-most-common "
                f"nt: {nt2ct}"
            )
            print(f"\t(Arbitrarily breaking tie: selecting alt = {alt_nt}.)")

    return cov, alt_freq, alt_nt


def get_alt_nt(pileup):
    """Raises an error if there are no mismatches.

    This returns the SECOND-MOST-COMMON nucleotide at a position, always.
    This second-most-common nucleotide will usually differ from the reference
    nucleotide at this position -- because we expect the reference nucleotide
    to match the consensus, right -- but this may not be the case. This is as
    expected: by always limiting this to the second-most-common nucleotide,
    we can ensure that freq(pos) = alt(pos) / reads(pos) remains within the
    easy-to-interpret range [0%, 50%]. For more detailed stuff on this, see the
    prokaryotic / eukaryotic classification paper...

    If you want to be able to exclude positions where the reference doesn't
    match the consensus, see get_alt_nt_if_reasonable().
    """
    _, alt_freq, alt_nt = get_alt_info_from_pleuk(pileup)

    if alt_freq == 0:
        raise ValueError("No mismatches at this position in the pileup.")

    return alt_nt
    

def is_reasonable(pileup):
    ri = pileup[1]
    return pileup[0][ri] == max(pileup[0])


def get_alt_nt_if_reasonable(pileup):
    """Like get_alt_nt(), but returns None in the case where the consensus
    (most-frequently-seen nucleotide at a position) does not match the
    reference. This also guarantees that the reference nucleotide will NOT be
    returned as the alt nucleotide, in the case of a tie.

    The rationale for this is that in some cases we want to treat a mismatch
    in a read as a "mutation", explicitly, from the reference. But if the
    reference doesn't match the consensus at this position, it becomes a bit
    weird trying to describe this as a mutation "from" something. I'm sure you
    could spend a lot of time arguing about the details of this case, but it's
    easiest to just exclude these rare positions.

    This won't return None if the consensus is tied for the
    most-frequently-seen nucleotide with another -- that case is ok. This
    matches the behavior of how we go through stuff for the codon mutation
    matrices.

    As with get_alt_nt(), this'll throw an error if there are no mismatches at
    this position.

    This code could probably be sped up a fair amount.
    """
    #_, alt_freq, alt_nt = get_alt_info_from_pleuk(pileup)
    cts = pileup[0]
    ri = pileup[1]
    if is_reasonable(pileup):
        # This position is "reasonable": the reference nucleotide is either
        # the most frequently seen nucleotide at this position, or tied for
        # this.
        nts = "ACGT"
        n2i = {"A": 0, "C": 1, "G": 2, "T": 3}
        # Slice out the reference nucleotide based on its index
        alt_nts = nts[:ri] + nts[ri + 1:]
        # Create alt_nt2ct, a dict mapping nucleotide to frequency at this
        # position -- this will only contain 3/4 of the nucleotides, with the
        # reference nucleotide excluded. Use this dict to figure out the
        # maximum-frequency alternate nucleotide (breaking ties arbitrarily)
        # and then return that.
        alt_nt2ct = {nt: cts[n2i[nt]] for nt in alt_nts}
        max_freq_alt_nt = max(alt_nt2ct, key=alt_nt2ct.get)
        # If all non-reference nucleotides have frequency 0, raise an error.
        if alt_nt2ct[max_freq_alt_nt] == 0:
            raise ValueError("No mismatches at this position in the pileup.")
        return max_freq_alt_nt
    else:
        # This position is "unreasonable", so just return None.
        return None
    

def get_mismatch_pcts(pileup):
    """Like get_mismatch_cts(), but returns percentages.

    The percentages include all aligned nucleotides (both matches and
    mismatches).
    """
    cov = get_cov(pileup, raise_error_if_0x=True)
    return [c / cov for c in get_mismatch_cts(pileup)]


def get_alt_nt_pct(pileup):
    """Returns the second-most-common nucleotide's relative freq at a position.

    This is the definition of alt(pos) used in the prok/euk classification
    report.

    So, e.g., if the pileup counts are [100, 20, 30, 50], then
    this'll return (50) / (100+20+30+50) = (50 / 200) = 0.25.

    If cov is 0, this'll just return 0. This behavior is debatable, but it
    should be sufficient for our purposes.
    """
    cov, alt_freq, alt_nt = get_alt_info_from_pleuk(pileup)
    if cov == 0:
        # Arguably this should be undefined, since it's n / 0 -- but
        # the desired result for these situations (we don't call a mutation
        # because this position is completely uncovered) is respected by
        # just treating the max-freq alternate nucleotide percentage as 0.
        return 0
    else:
        return alt_freq / cov


def naively_call_mutation(pileup, p, min_alt_pos=2, only_call_if_rare=False):
    """The main purpose of this file; naively "calls" a position as a
    p-mutation or not.

    A position pos with:
        second-most-common nt frequency = alt(pos),
        coverage = reads(pos)

    ...is a p-mutation
        (given some p in the range (0, 50] -- this corresponds to
        a percentage, so p = 50 indicates 50%)

    ...if (alt(pos) / reads(pos)) >= (p / 100).

    Or, as the paper puts it, freq(pos) >= p (replacing the left side of the
    equation with freq(pos), and expressing p as a value in (0, 0.5] instead).

    We previously had p passed to this function as a value in (0, 0.5], like in
    the paper, but in order to avoid the need to repeatedly divide p by 100
    before passing it to here -- and, simultaneously, to limit floating-point
    precision errors -- we now ask that p is passed in as a value in (0, 50].
    We handle things within this function, which makes life a bit easier.

    NOTE that this does not mean that p is now an integer, since you can
    totally have float percentages in this range: for example, 2.5 (indicating
    2.5%).

    Parameters
    ==========
    pileup: list
        A pileup entry for a given position. Same as the other pileup
        parameters throughout the functions in this file.

    p: float
        Should be in the range (0, 50]. Represents the threshold we use to
        call a mutation.

    min_alt_pos: int, default 2
        Should be a nonnegative integer (zero or one are ok, but I think
        they'll result in the same behavior). Even if a position qualifies as a
        p-mutation based on alt(pos) / reads(pos) >= p / 100, this position
        will NOT be considered a p-mutation unless a second condition,
        alt(pos) >= min_alt_pos, is met. The rationale for this is discussed in
        the paper.

    only_call_if_rare: bool, default False
        If True, then this will only call a p-mutation if
        (alt(pos) / reads(pos)) < HIGH_FREQUENCY_MIN_PCT / 100. The
        HIGH_FREQUENCY_MIN_PCT constant is defined above.

        This has the effect of allowing us to ignore "indisputable" mutations.

    Returns
    =======
    bool: True if this position is a p-mutation (and the min_alt_pos and
          only_call_if_rare conditions are met), False otherwise.
    """
    # Attempt to catch errors from me forgetting to update the definition of p
    # used throughout these analyses
    if p > 50 or p <= 0:
        raise ValueError(f"Hey p = {p} but it should be in the range (0, 50]")

    cov, alt_freq, alt_nt = get_alt_info_from_pleuk(pileup)
    return naively_call_mutation_directly(
        alt_freq, cov, p, min_alt_pos=min_alt_pos,
        only_call_if_rare=only_call_if_rare
    )


def naively_call_mutation_directly(
    alt_pos, cov_pos, p, min_alt_pos=2, only_call_if_rare=False
):
    """The guts of naively_call_mutation().

    Abstracted to a separate function so that I can call this from other
    contexts besides ordinary pileup -- e.g. for calling codons as p-mutated.
    """
    if p > 50 or p <= 0:
        raise ValueError(f"Hey p = {p} but it should be in the range (0, 50]")

    if alt_pos >= min_alt_pos:
        # We call a p-mutation if alt(pos) / reads(pos) >= p / 100.
        # Equivalently: we call a p-mutation if 100*alt(pos) >= p*reads(pos).
        #
        # This way, we avoid division; we aren't entirely out of the woods of
        # potential floating-point errors, but this should be more reliable, I
        # think. (And thanks to Python's support for arbitrary-precision
        # integers, we shouldn't need to worry about numbers getting
        # ridiculously large enough to cause overflow problems -- plus, like,
        # the the max value on the right hand side here would be what,
        # p=50 times a coverage of say 1,000,000x? That's no biggie.)
        lhs = 100 * alt_pos
        rhs = p * cov_pos
        if lhs >= rhs:
            # This position counts as a p-mutation, but we may still need to
            # make the "is this a rare mutation?" check.
            if only_call_if_rare:
                return is_position_rare_direct(alt_pos, cov_pos)
            else:
                return True
    else:
        return False


def is_position_rare(pileup):
    cov, alt_freq, alt_nt = get_alt_info_from_pleuk(pileup)
    return is_position_rare_direct(alt_freq, cov)


def is_position_rare_direct(alt_pos, cov_pos):
    lhs = 100 * alt_pos
    rhs_upper = HIGH_FREQUENCY_MIN_PCT * cov_pos
    return lhs < rhs_upper


def any_mismatches(pileup):
    """Returns True if alt(pos) > 0, False otherwise.

    Replicates the behavior of naive p-mutation calling with > p instead of
    the new >= p behavior at p = 0%.
    """
    _, alt_freq, _ = get_alt_info_from_pleuk(pileup)
    return alt_freq > 0


def get_deletions(pileup):
    return pileup[2]
