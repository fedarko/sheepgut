#! /usr/bin/env python3
# Given a BAM file listing alignments of reads to edges in a metagenome
# assembly graph, filters reads to just those almost fully mapped to a
# collection of "edges to focus on" (in our case, the three "genomes"
# discussed in our paper).
#
# This includes both filtering out "partially-mapped reads", which will
# mess with the mutation analyses, as well as just filtering out reads
# entirely(ish) mapped to other edges in the graph that we don't care about.
#
# THIS ASSUMES THAT SECONDARY ALIGNMENTS HAVE ALREADY BEEN FILTERED OUT.
# Supplementary alignments should have been left in, but secondary alignments
# will mess this up!!!
#
# This uses Pysam, in particular the fetch() function.
# bf.fetch(edge) returns an iterator over all AlignedSegments (alignments
# of [part of] reads) aligned to the sequence "edge" in the AlignmentFile "bf".
# This can include multiple reads with the same "query name," in the case
# of supplementary alignments being aligned to the same edge (e.g. a read
# from both the beginning and end of a genome).
# NOTE that a read may still have supplementary alignments to other edges
# besides the edge we're currently focusing on; we account for this in the
# pseudocode and code below.
#
# PSEUDOCODE:
# For each edge we want to focus on:
# | For each alignment of a read to this edge:
# | | Record the read's total length (if not already seen while looking for
# | | stuff in this edge).
# | |
# | | Record or update the number of (mis)matching positions in an alignment
# | | to the current edge that this read has; can be done by counting M/X/=
# | | operations in the CIGAR string (including both matches and mismatches).
# | | We can assume supplementary alignments to the same edge do not have
# | | overlapping coordinates, since we've already filtered the BAM to exclude
# | | these.
# |
# | For all other edges in this edge's weakly connected component, if
# | applicable:
# | | For each alignment of a read to this other edge:
# | | | If this read was aligned to the edge we're currently focusing on in
# | | | this component (i.e. we already saw above):
# | | | | Update the number of positions that this alignment of the read has
# | | | | by counting M/X/= operations, same as before.
# |
# | For each alignment of a read to the edge we want to focus on:
# | | Now that we've seen all alignments of this read in this component,
# | | compute the percentage of this read aligned to edges in this component.
# | | If this percentage passes a defined cutoff,
# | | write this alignment to an output BAM file; otherwise, don't. NOTE that
# | | this decision will be the same across all alignments of a given read
# | | (i.e. either all or none of the alignments of a read will be output to
# | | the BAM file).
#
# Then move on to the next edge to focus on and repeat the process, writing
# the kept reads to the same output BAM file as before. Notably, this could
# consider the same read more than once (i.e. if a read is aligned to
# multiple edges we intend to focus on) -- for two reasons. One, because
# specifically trying to avoid that case would get complicated, and two,
# because we could actually want that to happen if MIN_PERCENT_ALIGNED is
# small enough -- e.g. if we just want to keep reads 20% or more aligned to
# an edge, then the same read could totally be aligned 30% to one edge
# we want to focus on and 30% to another edge we want to focus on, etc.

import re
import pysam
from collections import defaultdict

# This is a percentage (value in [0, 1]). Reads where less than this
# percentage of the sequence is aligned to a given edge are filtered out of
# the BAM file.
MIN_PERCENT_ALIGNED = 0.9

print("Filtering out reads mapped to other edges and partially-mapped reads...")

edges_to_focus_on = ["edge_1671", "edge_2358", "edge_6104"]

# Edges 1671 and 2358 are in isolated components; however, edge 6104
# is in a component with 31 other edges (total 32 edges).
edges_in_ccs = {
    "edge_1671": ["edge_1671"],
    "edge_2358": ["edge_2358"],
    "edge_6104": [
        "edge_10465",
        "edge_11362",
        "edge_11363",
        "edge_11364",
        "edge_11365",
        "edge_11368",
        "edge_11369",
        "edge_11370",
        "edge_11371",
        "edge_11372",
        "edge_11375",
        "edge_11379",
        "edge_11380",
        "edge_11384",
        "edge_18000",
        "edge_18002",
        "edge_18003",
        "edge_18004",
        "edge_18005",
        "edge_18006",
        "edge_18007",
        "edge_18008",
        "edge_18009",
        "edge_18010",
        "edge_18087",
        "edge_18097",
        "edge_38487",
        "edge_6104",
        "edge_74421",
        "edge_9339",
        "edge_9340",
        "edge_9342"
    ]
}

# Input BAM file (contains read alignments to all edges in the graph)
bf = pysam.AlignmentFile("output/aln-sorted.bam", "rb")
# Output BAM file (will just contain alignments to the edges we're focusing on)
# These alignments will in turn be further filtered to alignments where a read
# was almost entirely mapped to a single edge or its component, to limit
# spurious mutations.
of = pysam.AlignmentFile("output/pmread-filtered-aln.bam", "wb", template=bf)

# Per the SAM v1 specification, top of page 8:
# "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ."
#
# These are the five CIGAR operations that consume character(s) from the query
# sequence (i.e. a read).
#
# We ignore S (soft clipping), since this indicates that a given position in a
# read is not matched to an edge; and we ignore I (insertion), since this also
# indicates that a position in a read is not really "matched" anywhere on the
# edge (although I guess you could argue that including insertions here may be
# useful; it probably isn't a big deal either way).
#
# By just looking at M, X, and = occurrences, we can get a count for each
# alignment of the number of bases "(mis)matched" to an edge.
# Since we have already filtered secondary alignments, and since per the SAM
# specification "A chimeric alignment is represented as a set of linear
# alignments that do not have large overlaps", we can sum these match counts
# across all alignments of a read and divide by the read length to get the
# approximate percentage of bases in the read aligned to an edge or group of
# edges.
#
# (Note that in practice it's possible for the supplementary alignments to
# share bases, in which case we could get percentages over 100%, but it's
# expected per the SAM spec that these overlaps should be small.)
matches = re.compile("(\d+)[MX=]")

def check_and_update_alignment(
    read, readname2len, readname2num_matches_in_cc, edge_name
):
    # Ensure that read length is consistent across all alignments involving
    # this read; also, make a record of previously unseen reads' lengths.
    # (Apparently query length is dependent on the actual alignment of this
    # read, so we use infer_read_length() instead because we care about the
    # actual length of the read)
    if read.query_name in readname2len:
        if readname2len[read.query_name] != read.infer_read_length():
            raise ValueError(
                "Inconsistent read lengths across alignments: {}: {}, {}".format(
                    read.query_name, 
                    readname2len[read.query_name],
                    read.infer_read_length()
                )
            )
    else:
        readname2len[read.query_name] = read.infer_read_length()

    # Each AlignedSegment returned by fetch(edge) should pertain to that
    # specific edge sequence -- this lets know that the CIGAR string of
    # this segment applies to the edge
    if read.reference_name != edge_name:
        raise ValueError(
            "Read reference name, {}, isn't {} as expected".format(
                read.reference_name, edge_name
            )
        )

    # The meat of this -- parse the CIGAR string of this alignment and
    # count all (mis)match operations, updating a defaultdict.
    allmatches = matches.findall(read.cigarstring)
    if allmatches:
        num_matches = sum([int(c) for c in allmatches])
        # Makes the simplifying assumption that supplementary alignments
        # do not overlap; otherwise we'd need to look at coordinates /
        # CIGAR operations and try to only count each position once
        readname2num_matches_in_cc[read.query_name] += num_matches
    else:
        # Raise an error if an alignment of this read does not involve
        # any (mis)match operations at all. This *could* happen in practice,
        # I guess, but if it does something is likely wrong. If this check
        # needs to be removed in the future, then this block could just be
        # commented out or replaced with a "pass" statement or something.
        raise ValueError(
            "No match chars (M/X/=) found in read {} CIGAR: {}".format(
                read.query_name, read.cigarstring
            )
        )

# Figure out all reads that are aligned to each edge to focus on
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
