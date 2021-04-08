# You can run these tests by just running "pytest"
from parse_pileup_read_bases import parse_read_bases


def test_one_indel():
    outs = parse_read_bases("act.acA.-2CGC,")
    assert len(outs) == 2
    # Three matches (. and ,)
    assert outs[0] == 3
    # Three unique non-matching nucleotides, ignoring case
    assert len(outs[1]) == 3
    # Of the non-matching bases, there are 3 A, 3 C, 1 T
    assert outs[1]["A"] == 3
    assert outs[1]["C"] == 3
    assert outs[1]["T"] == 1
    assert "G" not in outs[1]


def test_no_indels():
    outs = parse_read_bases("...,,ac...c")
    assert len(outs) == 2
    # 8 matches
    assert outs[0] == 8
    # two unique non-matching nucleotides, ignoring case
    assert len(outs[1]) == 2
    # Of the non-matching bases...
    assert outs[1]["A"] == 1
    assert outs[1]["C"] == 2
    assert "T" not in outs[1]
    assert "G" not in outs[1]


def test_many_indels():
    # "oops, all indels": everyone's fave breakfast cereal
    outs = parse_read_bases("+1AA-2CC.+3GTC")
    assert len(outs) == 2
    # 1 matches
    assert outs[0] == 1
    # one sneaky non-matching "A" after the +1A
    assert len(outs[1]) == 1
    assert outs[1]["A"] == 1


def test_only_indels():
    # "oops, all indels": everyone's fave breakfast cereal
    outs = parse_read_bases("+1A-2CC+3GTC")
    assert len(outs) == 2
    # no matches
    assert outs[0] == 0
    # no non-matching nucleotides
    assert len(outs[1]) == 0


def test_no_read_bases():
    # Works similarly to only indels case
    outs = parse_read_bases("")
    assert len(outs) == 2
    # no matches
    assert outs[0] == 0
    # no non-matching nucleotides
    assert len(outs[1]) == 0
