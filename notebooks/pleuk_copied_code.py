# NOTE: This file just contains a function copied from Pleuk's
# bam_utils library.
# I am duplicating this code here (rather than adding Pleuk as a dependency
# or submodule) because Pleuk is currently a private repo on GitHub, so it's
# annoying to integrate that with the SheepGut CI. It's easier for now to
# duplicate the code; when Pleuk is public, I can remove this duplicated code
# and instead just add a dependency on Pleuk.

def get_alt_pos_info(rec):
    """Returns info about the second-most-common nucleotide at a position.

    This nucleotide will usually differ from the reference nucleotide, but it
    may be the reference (i.e. at positions where the reference disagrees with
    the alignment's "consensus").

    This breaks ties arbitrarily.

    Parameters
    ==========
    rec: dict
        pysamstats record for a given position in an alignment produced
        by stat_variation().

    Returns
    =======
    (cov, alt nt freq, alt nt): tuple of (int, int, str)
        Describes the second-most-common nucleotide at a position.

        The first entry in this tuple is the (mis)match coverage at this
        position. This is an integer defined as the sum of A, C, G, T
        nucleotides at this position (note that this excludes degenerate
        nucleotides like N -- we could change this in the future if that'd be
        useful, I suppose). Note that this coverage could be zero, if no reads
        are aligned to this specific position.

        The second entry is the raw frequency of this nucleotide
        at this position: this will be an integer greater than or equal to 0.
        This is also referred to in the paper, etc. as alt(pos).

        The third entry is just the alternate nucleotide (one of A, C, G, T),
        represented as a string. This is returned for reference -- as of
        writing this isn't actually needed for Pleuk itself, but I have other
        code outside of Pleuk that benefits from this!
    """
    cov = rec["A"] + rec["C"] + rec["G"] + rec["T"]

    ordered_nts = sorted("ACGT", key=rec.get)

    # The literal nucleotide used in the numerator of freq(pos): one of A, C,
    # G, T
    alt_nt = ordered_nts[-2]

    # The raw frequency (in counts) of alt_nt. An integer >= 0.
    alt_nt_freq = rec[alt_nt]

    return (cov, alt_nt_freq, alt_nt)
