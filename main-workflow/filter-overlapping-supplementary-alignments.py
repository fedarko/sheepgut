#! /usr/bin/env python3
# Filters a BAM file so that reads with supplementary alignments that
# overlap with each other on the reference genome are removed.
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

import time
import pysam
from itertools import combinations
from collections import defaultdict

print("Filtering out overlapping supplementary alignments...")

t0 = time.time()

# Input BAM file (contains read alignments to all edges in the graph)
bf = pysam.AlignmentFile("output/aln-sorted.bam", "rb")


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


# Keeps track of the names of reads with overlapping supplementary alignments.
# We'll filter these reads out completely (so that they are not represented in
# any of the alignments remaining in the BAM file).
#
# I was initially going to maintain a sorted array for this, so that we could
# use binary search to quickly check if reads seen in the final "pass" over the
# BAM file had OSAs, but it turns out that in Python using sets is probably a
# better idea (or at the very least a simpler one):
# https://stackoverflow.com/a/212971
reads_with_osa = set()

for si, seq in enumerate(bf.references, 1):
    pct = 100 * (si / bf.nreferences)
    t1 = time.time()
    print(
        f"Pass 1/2: on seq {seq} ({si:,} / {bf.nreferences:,}) ({pct:.2f}%, "
        f"running for {t1 - t0:,.2f} sec...)"
    )

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

    # Identify overlapping alignments from the same read
    n_reads_w_osa_in_seq = 0
    for rn in readname2CoordsAndMQ:
        alns = readname2CoordsAndMQ[rn]
        if len(alns) > 1:
            # Okay, so this particular read has multiple supplementary
            # alignments to this sequence. Check if they overlap.

            for (a1, a2) in combinations(alns, 2):
                # Efficiently test for overlap between two ranges:
                # https://stackoverflow.com/a/3269471
                if a1[0] <= a2[1] and a2[0] <= a1[1]:
                    # Okay, these two alignments of this read overlap. 
                    reads_with_osa.add(rn)
                    n_reads_w_osa_in_seq += 1
                    break

    print(f"\t{n_reads_w_osa_in_seq:,} read(s) with overlapping supp alns...")

# Now we've made note of all reads with OSAs across *all* sequences in the
# alignments. We can make another pass through and output all reads without
# OSAs into a new BAM file.

# Output BAM file (filtered to remove reads with overlapping supplementary
# alignments, aka OSAs)
of = pysam.AlignmentFile(
    "output/overlap-supp-aln-filtered-aln.bam", "wb", template=bf
)

# TODO: maybe generalize this iteration code into a generator or something
# to limit code reuse
for si, seq in enumerate(bf.references, 1):
    pct = 100 * (si / bf.nreferences)
    t1 = time.time()
    print(
        f"Pass 2/2: on seq {seq} ({si:,} / {bf.nreferences:,}) ({pct:.2f}%, "
        f"running for {t1 - t0:,.2f} sec...)"
    )

    num_alns_retained = 0
    num_alns_filtered = 0
    for linearaln in bf.fetch(seq):
        rn = linearaln.query_name
        # If this read has OSAs anywhere in the alignment, don't include it in
        # the output BAM file. Otherwise, *do* include it!
        if rn in reads_with_osa:
            num_alns_filtered += 1
        else:
            of.write(linearaln)
            num_alns_retained += 1

    num_alns_total = num_alns_retained + num_alns_filtered
    if num_alns_total > 0:
        apct = 100 * (num_alns_retained / num_alns_total)
    else:
        apct = float("inf")
    print(
        f"\t{num_alns_retained:,} / {num_alns_total:,} ({apct:.2f}%) "
        "linear aln(s) retained."
    )

bf.close()
of.close()

t2 = time.time()
print("Filtered out overlapping supplementary alignments.")
print(f"Total time taken: {t2 - t0:,.2f} sec.")
