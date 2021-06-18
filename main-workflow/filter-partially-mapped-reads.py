#! /usr/bin/env python3
# Filters reads that are less than some percentage (MIN_PERCENT_ALIGNED)
# aligned to an edge in the graph, or to other edges adjacent to this edge in
# the repeat graph.
#
# Realized after the fact that this bears some resemblance to samclip
# (https://github.com/tseemann/samclip), although this differs a bit in the
# sort of alignments this allows to pass the filter and the sort of information
# it takes into account. These are probs ultimately minor distinctions, tho.

import re
import pysam
from collections import defaultdict
from utils import load_gfa

# This is a percentage (value in [0, 1]). Reads where less than this
# percentage of the sequence is aligned to a given edge are filtered out of
# the BAM file.
MIN_PERCENT_ALIGNED = 0.9

print("Filtering out reads mapped to other edges and partially-mapped reads...")

print("\tLoading graph...", end=" ")
with open("../config/input-graph", "r") as igfile:
    graph_filename = next(igfile).strip()

graph = load_gfa(graph_filename)
print("Done.")

# Input BAM file (contains alignments to all edges in the graph)
bf = pysam.AlignmentFile("output/overlap-supp-aln-filtered-and-sorted-aln.bam", "rb")
# Output BAM file (will just contain alignments that pass the filter)
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
# Since we have already filtered secondary alignments and overlapping
# supplementary alignments, we can sum these match counts across all
# alignments of a read and divide by the read length to get the approx
# percentage (see above for slight caveats) of bases in the read aligned to
# an edge or group of edges.
matches = re.compile("(\d+)[MX=]")

def check_and_update_alignment(
    aln, readname2len, readname2matchct, edge_name
):
    """Updates readname2len and readname2matchct based on an AlignedSegment.

    Note that readname2len and readname2matchct are both defined relative to
    a single edge in the graph (these aren't universal structures).
    """
    # Ensure that read length is consistent across all alignments involving
    # this read; also, make a record of previously unseen reads' lengths.
    # (Apparently query length is dependent on the actual alignment of this
    # read, so we use infer_read_length() instead because we care about the
    # actual length of the read)
    if aln.query_name in readname2len:
        if readname2len[aln.query_name] != aln.infer_read_length():
            raise ValueError(
                "Inconsistent read lengths across alignments: {}: {}, {}".format(
                    aln.query_name, 
                    readname2len[aln.query_name],
                    aln.infer_read_length()
                )
            )
    else:
        readname2len[aln.query_name] = aln.infer_read_length()

    # Each AlignedSegment returned by fetch(edge) should pertain to that
    # specific edge sequence -- this lets know that the CIGAR string of
    # this segment applies to the edge
    if aln.reference_name != edge_name:
        raise ValueError(
            "Alignment reference name, {}, isn't {} as expected".format(
                aln.reference_name, edge_name
            )
        )

    # The meat of this -- parse the CIGAR string of this alignment and
    # count all (mis)match operations, updating a defaultdict.
    allmatches = matches.findall(aln.cigarstring)
    if allmatches:
        matchct = sum([int(c) for c in allmatches])
        readname2matchct[aln.query_name] += matchct
    else:
        # Raise an error if an alignment of this read does not involve
        # any (mis)match operations at all. This *could* happen in practice,
        # I guess, but if it does something is likely wrong. If this check
        # needs to be removed in the future, then this block could just be
        # commented out or replaced with a "pass" statement or something.
        raise ValueError(
            "No match chars (M/X/=) found in read {} CIGAR: {}".format(
                aln.query_name, aln.cigarstring
            )
        )

# Figure out all reads that are aligned to each edge to focus on
for edge_to_focus_on in graph.nodes():
    print("Looking at edge {}...".format(edge_to_focus_on))
    # Maps read name to read length (which should be constant across all
    # alignments of that read). This variable is used both to store this info
    # (which is in turn used for sanity checking) as well as as a
    # crude indication of "have we seen this read yet?"
    readname2len = {}

    # Maps read name to number of match operations to the sequences of this
    # edge or adjacent edges in the graph.
    readname2matchct = defaultdict(int)

    for aln, i in enumerate(bf.fetch(edge_to_focus_on), 1):
        check_and_update_alignment(
            aln, readname2len, readname2matchct, edge_to_focus_on
        )
    print(f"\t{i} linear alignments to this edge.")
    
    # Identify adjacent edges to this one in the graph, if present. Allow
    # alignments to these edges to count towards readname2matchct (so we can
    # somewhat avoid penalizing edges that aren't in isolated components).
    # We could also expand this to include _all_ other edges in this edge's
    # component, but that will probs cause problems with hairball
    # component(s) in the graph that span thousands of edges...
    adj_edges = set(graph.neighbors(edge_to_focus_on)) - set([edge_to_focus_on])

    print(f"\t{len(adj_edges)} adjacent edges in the graph.")

    # Go through these "allowed" other edges; add to the number
    # of matches in cc for any reads that we see that we've already seen in
    # the edge to focus on. (We implicitly ignore any reads aligned to these
    # edges but not to the edge to focus on.)
    num_other_edge_alns_from_shared_reads = 0
    for other_edge in adj_edges:
        for aln in bf.fetch(other_edge):
            # If aln.query_name is NOT in readname2len, then this alignment's
            # read wasn't also aligned to the edge to focus on -- in this
            # case we implicitly ignore it, as mentioned above, since we only
            # really care about reads aligned to the edge we're focusing on.
            if aln.query_name in readname2len:
                check_and_update_alignment(
                    aln, readname2len, readname2matchct, other_edge
                )
                num_other_edge_alns_from_shared_reads += 1

    print(
        "f\t{num_other_edge_alns_from_shared_reads} linear alignments from "
        "shared reads to adjacent edges."
    )

    # Now that we've considered all relevant edges, we can compute
    # the approximate percentages of each read (aligned to the edge we're
    # focusing on) aligned to all edges in this component.
    #
    # And, finally, we can selectively write these reads' alignments
    # to an output BAM file accordingly.
    #
    # Note that although we're iterating over all alignments, these
    # computations are identical for alignments from the same read (at least in
    # the context of the same edge_to_focus_on). So if a read passes the
    # filter for an edge, all its alignments to this edge will be output to the
    # BAM file; and if the read fails the filter for an edge, none of its
    # alignments to that edge will be output to the BAM file (although this
    # does not preclude other alignments from this read to other edges from
    # being output in the context of other edges in this script).
    p = 0
    for aln in bf.fetch(edge_to_focus_on):
        if aln.query_name not in readname2len:
            raise ValueError(
                f"We should have seen read {aln.query_name} earlier!"
            )
        read_len = readname2len[aln.query_name]

        if readname2matchct[aln.query_name] == 0:
            # As with a similar error case in check_and_update_alignment(),
            # I *guess* this could happen in practice but it is probably
            # indicative of an error in most cases.
            raise ValueError(
                f"Read {aln.query_name} had no match operations done?"
            )

        read_matchct_in_cc = readname2matchct[aln.query_name]
        perc = read_matchct_in_cc / read_len

        if perc >= MIN_PERCENT_ALIGNED:
            of.write(read)
            p += 1
    print(
        f"\t{p} / {i} ({((p / i)*100):.2f}%) of alignments in "
        f"edge {edge_to_focus_on} passed the filter."
    )

bf.close()
of.close()

print("Filtered to fully (ish) aligned reads.")
