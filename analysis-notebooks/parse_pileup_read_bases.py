import re
from collections import defaultdict


def parse_read_bases(read_bases_column):
    """Given a str defining a "read bases" column in a pileup file, returns a
    2-tuple of (matchct, non_matches).

    matchct is an int indicating the number of matches present at this
    position.

    non_matches is a dict mapping nucleotides that do not match the
    reference base (case-insensitive) to the number of occurrences present
    at this position.
    """
    # Remove all "^X" (where X is any character) matches in the read
    # bases. The reason for doing this is that the X indicates mapping
    # quality (see http://samtools.sourceforge.net/pileup.shtml), and X
    # can be a "." or a "," -- which will break our code when we try to
    # count matches (which are represented in the read bases as "." or
    # ","). Failing to do this step will result in, rarely, some lines
    # having more matches than coverage, due to these "^X" patterns.
    rb_no_mq = re.sub(r"\^.", "", read_bases_column)

    matchct = 0
    non_matches = defaultdict(int)

    chars = list(rb_no_mq)
    i = 0
    while i < len(chars):
        if chars[i] == "-" or chars[i] == "+":
            # This is an indel -- see how long it is, and skip forward that
            # many characters in the column. This approach roughly based on
            # https://stackoverflow.com/a/40231190, insofar as we first parse
            # the indel length and then use that to determine how far forward
            # to go.
            indel_len_match = re.match(r"(\d+)", rb_no_mq[i + 1 :])
            if indel_len_match is None:
                # We just found a - or + that wasn't followed by at least one
                # number. This indicates that this is a malformed read bases
                # column (we already filtered out the only ASCII quality stuff
                # that should be present in this column).
                raise ValueError("- or + followed by no number")
            else:
                indel_len_str = indel_len_match.group(1)
                indel_len = int(indel_len_str)
                # Jump forward by 1 + len(indel_len_str) + indel_len steps.
                #
                #  1: We add an extra +1 to skip over the - or + that starts
                #  the "declaration" of this indel.
                #
                #  len(indel_len_str): We skip over the characters that
                #  describe the length of the indel.
                #
                #  indel_len: We skip over the actual characters in the indel.
                #
                # For example, we may run into something like "..+2CGC,..".
                # When i = 2 (i.e. we're at the "+"), we should jump forward by
                # 1 + 1 + 2 = 4 characters: 1 character to get past the "+", 1
                # character to get past the "2", and 2 characters to get past
                # the "CG". So now i = 6, and we're right where we need to be
                # (at the "C").
                i += 1 + len(indel_len_str) + indel_len
                continue

        elif chars[i] == "." or chars[i] == ",":
            matchct += 1

        elif chars[i] in "ACGTacgt":
            non_matches[chars[i].upper()] += 1

        elif chars[i] in "Nn":
            # This should never happen unless we're working with
            # sequences/reads with gaps. Which I'm not.
            raise ValueError(
                "Lemme get uhhh gapless assembly with uhhhhhhh 2 liter "
                "of Campylobacter"
            )

        # Jump forward 1 character (going to the next thing in the read bases)
        i += 1

    return (matchct, non_matches)
