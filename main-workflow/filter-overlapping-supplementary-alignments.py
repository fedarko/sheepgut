#! /usr/bin/env python3
# Filters a BAM file so that no supplementary alignments of a given read
# overlap with other supplementary alignments of the same read.
#
# This is applied to the entire BAM file (and since this is done before
# the partially-mapped-read filtering script, this therefore is applied
# to all alignments in the BAM file -- not just to the edges of interest,
# or to other edges in their components). Could be sped up, of course, to only
# be applied to edges in these components.
#
# THIS ASSUMES THAT SECONDARY ALIGNMENTS HAVE ALREADY BEEN FILTERED OUT.
# Supplementary alignments should have been left in, but secondary alignments
# will mess this up!!!

import re
import pysam
import networkx as nx
from itertools import combinations
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

def get_quartet(alnseg):
    """Returns a tuple of information distinguishing a linear alignment.

    Parameters
    ----------
        alnseg: pysam.AlignedSegment
    
    Returns
    -------
        (s, e, mq, st): (int, int, int, str)
            Segment start, end, mapping quality, and to_string() output.

            The start and end are both inclusive, to simplify comparison of
            alignment ranges for detecting overlaps. 

            The reason we include the fourth element (to_string()) is to make
            it easier to distinguish linear alignments from the same read. It
            is very unlikely (but still possible I guess) that multiple
            alignments from a read will have identical QUAL values AND
            identical tags, both of which are included in to_string() as of
            writing.

    Raises
    ------
        ValueError
            If the segment's start is greater than its end (both in inclusive
            coordinates). (If this ends up being a problem in practice, maybe
            because there of reverse-mapped reads or something (???), then this
            could probs be modified to just reverse the start and end in this
            case.)
    """
    # We add 1 to the end since this is a half-open interval -- we want
    # the coordinates we use for computing overlap to be completely
    # inclusive intervals
    s = linearaln.reference_start
    e = linearaln.reference_end + 1
    if s > e:
        raise ValueError(
            f"Malformed linear alignment coordinates: start {s}, end {e}"
        )
    mq = linearaln.mapping_quality
    st = linearaln.to_string()
    return (s, e, mq, st)

for si, seq in enumerate(bf.references, 1):
    pct = 100 * (si / bf.nreferences)
    print(f"On seq {seq} ({si:,} / {bf.nreferences:,}) ({pct:.2f}%)...")

    # Identify all linear alignments of each read to this sequence
    readname2CoordsAndMQ = defaultdict(list)
    num_linear_alns = 0
    for linearaln in bf.fetch(seq):
        rn = linearaln.query_name
        alndetails = get_quartet(linearaln)
        if alndetails in readname2CoordsAndMQ[rn]:
            raise ValueError(
                f"Indistinguishable linear alignments to seq {seq} with read "
                f"name {rn}: multiple reads share (start, end, mapq, "
                f"to_string()) of {alndetails}"
            )
        readname2CoordsAndMQ[rn].append(alndetails)
        num_linear_alns += 1

    # The number of unique reads is just the number of keys in this dict
    num_reads = len(readname2CoordsAndMQ)

    print(f"\t{num_reads:,} unique read(s), {num_linear_alns:,} linear aln(s)...")

    # Identify and remove overlapping linear alignments from the same read
    num_reads_with_osa = 0
    for rn in readname2CoordsAndMQ:
        alns = readname2CoordsAndMQ[rn]
        if len(alns) > 1:
            # Okay, so this particular read has multiple supplementary
            # alignments to this sequence. Check if they overlap.
            # We model this by constructing a graph where nodes correspond to
            # linear alignments of this read, and edges connect overlapping
            # alignments.
            ag = nx.Graph()

            # add nodes
            for a in alns:
                ag.add_node(a)

            # add edges
            for (a1, a2) in combinations(alns, 2):
                # Efficiently test for overlap between two ranges:
                # https://stackoverflow.com/a/3269471
                if a1[0] <= a2[1] and a2[0] <= a1[1]:
                    # Okay, these two alignments of this read overlap. 
                    ag.add_edge(a1, a2)
            
            if len(ag.edges) > 0:
                num_reads_with_osa += 1

            # Remove alignments until no overlaps remain.
            #
            # There are many possible ways to do this, depending on what we
            # prioritize. For example, we may want to see if there are any
            # alignments with many edges in the graph, and remove these
            # alignments first -- to limit the total amount of alignments we
            # need to remove (imagine the graph of a1 -- a2 -- a3; we could
            # remove a1 and a3, or we could remove just a2).
            #
            # However, here we take a simpler approach and just assume that
            # most nodes will have degree 1. So we consider each edge in
            # isolation and just remove the alignment on this edge
            # with the lower mapping quality.
            while len(ag.edges) > 0:
                arbitrary_edge = list(ag.edges)[0]
                # compare alignment mapping qualities: higher is better.
                # see https://samtools.github.io/hts-specs/SAMv1.pdf
                if arbitrary_edge[0][2] > arbitrary_edge[1][2]:
                    to_remove = arbitrary_edge[0]
                else:
                    # note that this also includes the case where the mapping
                    # qualities are equal (in that case, the decision is really
                    # arbitrary, I guess? we could do things like consider
                    # alignment *spans* [i.e. end - start + 1 or something] but
                    # that could get misleading if some alignments have a lot
                    # of skips, etc. Could also consider number of matches in
                    # CIGAR string, too...? If we wanted to get fancy.)
                    to_remove = arbitrary_edge[1]

                ag.remove_node(to_remove)
                readname2CoordsAndMQ[rn].remove(to_remove)
    print(f"\t{num_reads_with_osa:,} read(s) with overlapping supp alns...")

    # Now, go and write out only the linear alignments we retained
    num_alns_retained = 0
    for linearaln in bf.fetch(seq):
        rn = linearaln.query_name
        alndetails = get_quartet(linearaln)
        if alndetails in readname2CoordsAndMQ[rn]:
            of.write(linearaln)
            num_alns_retained += 1
    print(f"\t{num_alns_retained:,} linear aln(s) retained.")

bf.close()
of.close()

print("Filtered out overlapping supplementary alignments.")
