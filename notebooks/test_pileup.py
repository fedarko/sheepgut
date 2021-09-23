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
