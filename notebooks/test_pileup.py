import pytest
import pileup

def test_get_cov():
    assert pileup.get_cov([[10, 0, 3, 50], 3, 0]) == 63
    assert pileup.get_cov([[10, 0, 3, 50], 0, 0]) == 63
    assert pileup.get_cov([[0, 0, 0, 0], 1, 0]) == 0
    with pytest.raises(ValueError, match="has coverage of 0x"):
        pileup.get_cov([[0, 0, 0, 0], 1, 0], raise_error_if_0x=True)


def test_get_alt_info_from_pleuk():
    # The reference SHOULD NOT make a difference in how the alt nucleotide is
    # determined
    for ri in (0, 1, 2, 3):
        assert pileup.get_alt_info_from_pleuk([[10, 0, 3, 50], ri, 0]) == (63, 10, "A")

    ai = pileup.get_alt_info_from_pleuk([[10, 10, 10, 10], 0, 0])
    assert ai[0] == 40
    assert ai[1] == 10
    # ties are broken arbitrarily
    assert ai[2] in "ACGT"

    ai2 = pileup.get_alt_info_from_pleuk([[5, 10, 10, 0], 0, 0])
    assert ai2[0] == 25
    assert ai2[1] == 10
    # ties are broken arbitrarily
    assert ai2[2] in "CG"


def test_get_alt_nt_pct():
    # 50 / 200 = 0.25
    assert pileup.get_alt_nt_pct([[100, 20, 30, 50], 0, 5]) == 0.25

    # silly case: all zeroes
    assert pileup.get_alt_nt_pct([[0, 0, 0, 0], 0, 5]) == 0


def test_get_alt_nt():
    # This doesn't care about the reference nucleotide at this position, so --
    # if the reference is NOT the consensus -- this can return the reference as
    # the "alt" nt. This is as expected.
    for ri in (0, 1, 2, 3):
        assert pileup.get_alt_nt([[100, 20, 30, 50], ri, 5]) == "T"

    # Break ties arbitrarily
    for ri in (0, 1, 2, 3):
        assert pileup.get_alt_nt([[100, 200, 200, 50], ri, 5]) in "CG"
        assert pileup.get_alt_nt([[200, 200, 300, 50], ri, 5]) in "AC"
        assert pileup.get_alt_nt([[100, 100, 100, 100], ri, 5]) in "ACGT"

    # Raise an error in the case where no alt nucleotides are present.
    # This error also applies to 0x coverage positions!
    bad_cases = (
        [[5, 0, 0, 0], 0, 3],
        [[0, 0, 0, 0], 0, 3]
    )
    for pu in bad_cases:
        with pytest.raises(
            ValueError, match=r"No mismatches at this position in the pileup."
        ):
           pileup.get_alt_nt(pu)


def test_get_alt_nt_if_reasonable():
    # Unlike get_alt_nt(), this does take the reference nt into account. If the
    # reference is the consensus, this behaves as expected; if not, then
    # this'll just return None.
    assert pileup.get_alt_nt_if_reasonable([[100, 20, 30, 50], 0, 5]) == "T"
    for ri in (1, 2, 3):
        assert pileup.get_alt_nt_if_reasonable([[100, 20, 30, 50], ri, 5]) is None

    # Even if the reference nt is tied for the most frequent nt, it should
    # never be returned as the alt nt by this function.
    assert pileup.get_alt_nt_if_reasonable([[100, 200, 200, 50], 0, 5]) is None
    assert pileup.get_alt_nt_if_reasonable([[100, 200, 200, 50], 1, 5]) == "G"
    assert pileup.get_alt_nt_if_reasonable([[100, 200, 200, 50], 2, 5]) == "C"
    assert pileup.get_alt_nt_if_reasonable([[100, 200, 200, 50], 3, 5]) is None

    # Ties are still broken arbitrarily if there are multiple most-frequent
    # non-reference nucleotides
    assert pileup.get_alt_nt_if_reasonable([[200, 200, 300, 50], 0, 5]) is None
    assert pileup.get_alt_nt_if_reasonable([[200, 200, 300, 50], 1, 5]) is None
    assert pileup.get_alt_nt_if_reasonable([[200, 200, 300, 50], 2, 5]) in "AC"
    assert pileup.get_alt_nt_if_reasonable([[200, 200, 300, 50], 3, 5]) is None

    assert pileup.get_alt_nt_if_reasonable([[100, 100, 100, 100], 0, 5]) in "CGT"
    assert pileup.get_alt_nt_if_reasonable([[100, 100, 100, 100], 1, 5]) in "AGT"
    assert pileup.get_alt_nt_if_reasonable([[100, 100, 100, 100], 2, 5]) in "ACT"
    assert pileup.get_alt_nt_if_reasonable([[100, 100, 100, 100], 3, 5]) in "ACG"

    # Same error-checking as in get_alt_nt()
    bad_cases = (
        [[5, 0, 0, 0], 0, 3],
        [[0, 0, 0, 0], 0, 3]
    )
    for pu in bad_cases:
        with pytest.raises(
            ValueError, match=r"No mismatches at this position in the pileup."
        ):
           pileup.get_alt_nt_if_reasonable(pu)


def test_naively_call_mutation():
    # Mutation calling, since we're computing alt(pos) differently now (based
    # on 2nd-most-common nt, rather than most-common non-reference nt), is not
    # impacted by the reference nt. We test all possible reference nt values to
    # verify that the results are consistent.
    for ri in (0, 1, 2, 3):

        # TEST CASE 1: 100 G at this position; 0 of all other nts
        not_mutated = [[0, 0, 100, 0], ri, 0]
        assert not pileup.naively_call_mutation(not_mutated, 0)

        # TEST CASE 2: 97 G at this position; 1 C; 0 A; 0 T
        mutated = [[0, 3, 97, 0], ri, 0]
        # test values of p for which this is a p-mutation
        for mp in (0, 0.01, 0.02, 0.025):
            assert pileup.naively_call_mutation(mutated, mp)
        # test values of p for which this is not a p-mutation
        for nmp in (0.03, 0.035, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
            assert not pileup.naively_call_mutation(mutated, nmp)

    tv = [[0, 0, 100, 0], 2, 0]
    for bad_p in (0.51, 0.6, 0.7, 0.8, 0.9, 1, 50, 100, -0.01, -10, -100):
        with pytest.raises(
            ValueError, match=r"should be in the range \[0, 0.5\]"
        ):
            pileup.naively_call_mutation(tv, bad_p)
