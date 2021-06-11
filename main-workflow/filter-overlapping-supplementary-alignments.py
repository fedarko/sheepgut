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

def get_triplet(alnseg):
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
    return (s, e, mq)

for si, seq in enumerate(bf.references, 1):
    print(f"On seq {seq} ({si:,} / {bf.nreferences:,})...")

    # Identify all linear alignments of each read to this sequence
    readname2CoordsAndMQ = defaultdict(list)
    num_linear_alns = 0
    for linearaln in bf.fetch(seq):
        rn = linearaln.query_name
        alndetails = get_triplet(linearaln)
        if alndetails in readname2CoordsAndMQ[rn]:
            raise ValueError(
                f"Indistinguishable linear alignments to seq {seq} with read "
                f"name {rn}: multiple reads share (start, end, mapq) of "
                f"{alndetails}"
            )
        readname2CoordsAndMQ[rn].append(alndetails)
        num_linear_alns += 1

    # The number of unique reads is just the number of keys in this dict
    num_reads = len(readname2CoordsAndMQ)

    print(f"\t{num_reads:,} unique reads, {num_linear_alns:,} linear alns...")

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
                    # of skips, etc.)
                    to_remove = arbitrary_edge[1]

                ag.remove_node(to_remove)
                readname2CoordsAndMQ[rn].remove(to_remove)
    print(f"\t{num_reads_with_osa:,} reads with overlapping supp alns...")

    # Now, go and write out only the linear alignments we retained
    num_alns_retained = 0
    for linearaln in bf.fetch(seq):
        rn = linearaln.query_name
        alndetails = get_triplet(linearaln)
        if alndetails in readname2CoordsAndMQ[rn]:
            of.write(linearaln)
            num_alns_retained += 1
    print(f"\t{num_alns_retained:,} linear alignments retained.")

bf.close()
of.close()

print("Filtered out overlapping supplementary alignments.")

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
