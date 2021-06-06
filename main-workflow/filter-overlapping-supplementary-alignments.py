#! /usr/bin/env python3
# Filters a BAM file so that no supplementary alignments of a given read
# overlap with other supplementary alignments of the same read.
#
# This is applied to the entire BAM file (and since this is done before
# the partially-mapped-read filtering script, this therefore is applied
# to all alignments in the BAM file -- not just to the edges of interest,
# or to other edges in their components). Could be sped up, of course.
#
# THIS ASSUMES THAT SECONDARY ALIGNMENTS HAVE ALREADY BEEN FILTERED OUT.
# Supplementary alignments should have been left in, but secondary alignments
# will mess this up!!!

import re
import pysam
import networkx
from collections import defaultdict

print("Filtering out overlapping supplementary alignments...")

# Input BAM file (contains read alignments to all edges in the graph)
bf = pysam.AlignmentFile("output/aln-sorted.bam", "rb")

# Output BAM file (will just contain alignments to the edges we're focusing on)
# These alignments will in turn be further filtered to alignments where a read
# was almost entirely mapped to a single edge or its component, to limit
# spurious mutations.
of = pysam.AlignmentFile(
    "output/overlap-supp-aln-filtered-aln.bam", "wb", template=bf
)

# -For all seqs in the BAM file
#  -Set up defaultdict(list): readname2coordsAndMapQ
#  -Set up defaultdict(int): readname2numSeqAlns
#  -For all linear alignments of reads to this seq
#   -Record the coordinates of this alignment, and its mapping quality, in the
#    dict. Should be a list of 3-tuples, e.g. (0, 100, 240) for a read that
#    goes from 0 to 99 with MAPQ of 240.
#  -Consider all reads where readname2numSeqAlns > 1
#   -See if any of these reads' alignments overlap with each other.
#    "Overlap" is going to be sort of an all vs. all relationship btwn. all
#    alignments of this read to this sequence -- think of it like a graph,
#    let's say with 0-indexed locations in readname2coordsAndMapQ indicating
#    nodes for linear alignments. An edge btwn. two nodes indicates overlap.
#    Note that if 0 and 1 overlap, and 1 and 2 overlap, this doesn't
#    necessarily mean that 0 and 2 overlap.
#   -So we construct this graph, and then go through all edges in it. We want
#    to minimize the number of supplementary alignments we have to throw out
#    while maximizing the total MAPQ of all the alignments we leave in.
#    There is probably a fancy way to solve this, but we can use a simple
#    greedy algorithm where we consider each edge and remove the alignment with
#    the worst mapping quality on the edge (choosing arbitrarily if both
#    alignments on the edge have equal MAPQ). When all edges are removed, we
#    are done.
#   -We then write to the output BAM file all of the alignments we DIDN'T
#    remove.

for edge_to_focus_on in edges_to_focus_on:
    print("Looking at edge {}.".format(edge_to_focus_on))
    # Maps read name to read length (which should be constant across all
    # alignments of that read). This variable is used both to store this info
    # and as a crude indication of "have we seen this read yet?"
    readname2len = {}

    # Maps read name to number of match operations to the sequences of edges
    # in edges_in_ccs[edge_to_focus_on].
    readname2num_matches_in_cc = defaultdict(int)

    i = 0
    for read in bf.fetch(edge_to_focus_on):
        check_and_update_alignment(
            read, readname2len, readname2num_matches_in_cc, edge_to_focus_on
        )
        i += 1
    print("{} alignments for this edge.".format(i))
    # Go through other edges in this component, if present; add to the number
    # of matches in cc for any reads that we see that we've already seen in
    # the edge to focus on. (We implicitly ignore any reads aligned to these
    # edges but not to the edge to focus on.)
    other_edges = set(edges_in_ccs[edge_to_focus_on]) - set([edge_to_focus_on])
    for other_edge in other_edges:
        for read in bf.fetch(other_edge):
            # If read.query_name is NOT in readname2len, then this read wasn't
            # aligned to the edge to focus on -- in this case we implicitly
            # ignore it, as mentioned above, since we only really care about
            # reads aligned to the edge we're focusing on.
            if read.query_name in readname2len:
                check_and_update_alignment(
                    read, readname2len, readname2num_matches_in_cc, other_edge
                )

    # Now that we've considered all edges in this component, we can compute
    # the approximate percentages of each read (aligned to the edge we're
    # focusing on) aligned to all edges in this component.
    #
    # And, finally, we can selectively write these reads to an output BAM file
    # accordingly.
    #
    # NOTE that "read" is kind of a misleading variable name here, since
    # the same read can be listed multiple times in a BAM/SAM file -- we're
    # really iterating over alignments, where the same read could be aligned
    # multiple times. In this case we either output all of these alignments for
    # a given read (if that read passes the percentage filter) or none of them
    # (if that read does not pass the percentage filter).
    p = 0
    n = 0
    for read in bf.fetch(edge_to_focus_on):
        if read.query_name not in readname2len:
            raise ValueError("We should have seen this read earlier!")
        read_len = readname2len[read.query_name]

        if readname2num_matches_in_cc[read.query_name] == 0:
            # As with a similar error case in check_and_update_alignment(),
            # I *guess* this could happen in practice but it is probably
            # indicative of an error in most cases.
            raise ValueError("This read had no match operations done?")

        read_num_matches_in_cc = readname2num_matches_in_cc[read.query_name]

        perc = read_num_matches_in_cc / read_len
        # small sanity check, print out first 10 alignments for each read
        if n < 10:
            print(
                "FYI: alignment of read {} has {} matches, len {}, {}%".format(
                    read.query_name, read_num_matches_in_cc, read_len,
                    perc * 100
            ), end="")
        if perc >= MIN_PERCENT_ALIGNED:
            of.write(read)
            if n < 10:
                print("; passed!")
            p += 1
        else:
            if n < 10:
                print("; failed!")
        n += 1
    print(
        "{} / {} ({:.2f}%) of alignments in edge {} passed the filter.".format(
            p, i, (p/i) * 100, edge_to_focus_on
        )
    )

bf.close()
of.close()

print("Filtered to fully (ish) aligned reads.")
