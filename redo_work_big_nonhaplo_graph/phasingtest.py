#! /usr/bin/env python3
import statistics
import pysam
import skbio
from itertools import combinations

bf = pysam.AlignmentFile("aln-sorted.bam")

# The gene's coordinates (1-indexed) are [1,208,927, 1,210,075] (that's an
# inclusive range).
# For 0-indexing in Python / Pysam notation, we represent this as
# [1,208,926, 1,210,075) -- Python / Pysam end coordinates are not included,
# so this really ranges to 1,210,074.
g1 = (1208926, 1210075)

# Figure out all reads that completely cover this gene (i.e. they are aligned
# to start before/at the start of the gene and end after/at the end of the
# gene).
covering_reads = []
for read in bf.fetch("edge_6104", g1[0], g1[1]):
    # fetch() returns all reads that are incident on a region, but this
    # includes reads that don't fully cover the region. Hence our checking that
    # the read covers the gene on both sides (it doesn't start and/or end
    # within the middle of the gene).
    read_pos = read.get_reference_positions()
    if read_pos[0] <= g1[0] and read_pos[-1] >= g1[1] - 1:
        covering_reads.append(read)

# These are the 34 notably mutated positions in this gene. They're 1-indexed.
mutated_positions_1indexed = [1209001, 1209010, 1209022, 1209058, 1209104, 1209115, 1209121, 1209126, 1209133, 1209136, 1209142, 1209145, 1209148, 1209154, 1209205, 1209241, 1209265, 1209266, 1209297, 1209322, 1209325, 1209337, 1209400, 1209403, 1209418, 1209421, 1209424, 1209523, 1209577, 1209610, 1209697, 1209700, 1209757, 1209796]
# Convert to 0-indexing.
mutated_positions = [p - 1 for p in mutated_positions_1indexed]

# Extract the reference sequence of the gene. This will let us figure out
# easily (ish) which reads are mutated at which positions.
camp = skbio.DNA.read("edge_6104.fasta")
geneseq = camp[g1[0]:g1[1]]

# Go through all pairs of mutated positions (ignoring ordering). There will be
# (34 choose 2) = 561 total pairs.
num_bothmutated_reads_across_all_mutated_pos_pairs = []
ii = 1
for mpospair in combinations(mutated_positions, 2):
    # Save the mutated positions in easier-to-handle variables.
    p0, p1 = mpospair
    if p0 == p1: raise ValueError("p0 should never equal p1")
    if p0 > p1: raise ValueError("p0 should be < p1")
    # Find out the reference nucleotides at p0 and at p1.
    p0ref = str(camp[p0])
    p1ref = str(camp[p1])
    print("Pair {} / 561. Ref p0 = {}, p1 = {}; p0 = {}, p1 = {}".format(
        ii, p0, p1, p0ref, p1ref
    ))
    # Figure out how many covering reads have mutations at BOTH p0 and p1.
    num_bothmutated_reads = 0
    for read in covering_reads:
        # Figure out where, exactly, in this read was aligned to p0 and p1.
        ap = read.get_aligned_pairs()
        query_pos_aligned_to_p0 = None
        query_pos_aligned_to_p1 = None
        for pair in ap:
            # (... this assumes that p0 will never equal p1. that'd be bad.)
            if pair[1] == p0:
                query_pos_aligned_to_p0 = pair[0]
            elif pair[1] == p1:
                query_pos_aligned_to_p1 = pair[0]
                # assumes p0 < p1; lets us avoid unnecessary iteration
                break
        # If this read wasn't aligned to either of these positions (e.g. it was
        # skipped), or if there was an insertion / deletion in the alignment,
        # then just ignore this read.
        if query_pos_aligned_to_p0 is None or query_pos_aligned_to_p1 is None:
            continue
        # Now that we know the exact positions in this read that were aligned
        # to p0 and p1, extract the nucleotides at these positions in the read
        p0query = read.query_sequence[query_pos_aligned_to_p0]
        p1query = read.query_sequence[query_pos_aligned_to_p1]
        # Does the read disagree at BOTH positions?
        if p0query != p0ref and p1query != p1ref:
            #print("Disagreement! Read name = {}".format(read.query_name))
            #print("Read p0 pos = {}, seq = {}".format(query_pos_aligned_to_p0, p0query))
            #print("Read p1 pos = {}, seq = {}".format(query_pos_aligned_to_p1, p1query))
            num_bothmutated_reads += 1
    num_bothmutated_reads_across_all_mutated_pos_pairs.append(
        num_bothmutated_reads
    )
    print("Pair of mutations {} had {} reads mutated at both p0 and p1.".format(
        mpospair, num_bothmutated_reads
    ))
    ii += 1

print("Average # of both-mutated reads across all mutated position pairs: {}".format(statistics.mean(num_bothmutated_reads_across_all_mutated_pos_pairs)))
print("Std. dev. of both-mutated reads across all mutated position pairs: {}".format(statistics.stdev(num_bothmutated_reads_across_all_mutated_pos_pairs)))
