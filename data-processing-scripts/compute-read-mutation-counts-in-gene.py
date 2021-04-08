#! /usr/bin/env python3
import pysam
import skbio

bf = pysam.AlignmentFile("aln-sorted.bam")

# The gene's coordinates (1-indexed) are [1,208,927, 1,210,075] (that's an
# inclusive range).
# For 0-indexing in Python / Pysam notation, we represent this as
# [1,208,926, 1,210,075) -- Python / Pysam end coordinates are not included,
# so this really ranges to 1,210,074.
g1 = (1208926, 1210075)

camp = skbio.DNA.read("edge_6104.fasta")

# Figure out all reads that completely cover this gene (i.e. they are aligned
# to start before/at the start of the gene and end after/at the end of the
# gene).
num_covering_reads = 0
num_mutations_in_gene_per_read = []
for read in bf.fetch("edge_6104", g1[0], g1[1]):
    # fetch() returns all reads that are incident on a region, but this
    # includes reads that don't fully cover the region. Hence our checking that
    # the read covers the gene on both sides (it doesn't start and/or end
    # within the middle of the gene).
    read_pos = read.get_reference_positions()
    if read_pos[0] <= g1[0] and read_pos[-1] >= g1[1] - 1:
        num_mutations_in_gene = 0
        num_covering_reads += 1
        ap = read.get_aligned_pairs()
        for pair in ap:
            if pair[1] in range(g1[0], g1[1]) and pair[0] is not None:
                query_pos = pair[0]
                query_seq = read.query_sequence[query_pos]
                ref_seq = str(camp[pair[1]])
                if query_seq != ref_seq:
                    num_mutations_in_gene += 1
        num_mutations_in_gene_per_read.append(num_mutations_in_gene)
        if num_covering_reads > 500:
            print("{} covering reads so far".format(num_covering_reads))

print("{} covering reads.".format(num_covering_reads))
print("{}".format(sorted(num_mutations_in_gene_per_read)))
