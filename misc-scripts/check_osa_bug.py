#! /usr/bin/env python3
# Determine how many reads were erroneously filtered by the off-by-2 error in
# the OSA filter. This code was scavenged from strainFlye's align_utils code.

import pysam
from itertools import combinations
from collections import defaultdict
from strainflye import cli_utils, bam_utils
from strainflye.errors import WeirdError


def check_osa_bug(in_bam, out_reads_info):
    """Identifies reads that were incorrectly flagged as having OSAs.

    This happened if a read had no OSAs ordinarily, but if any two alignments
    of the read were within exactly one or two nts away from each other. The
    problem was that, rather than subtracting 1 from the end coord for each
    linear alignment, I was adding 1, causing an off-by-2 error.

    Parameters
    ----------
    in_bam: str
        Location of the BAM file to check.

    out_reads_info: str
        File to which info about impacted reads will be written.
    """
    fancylog = cli_utils.fancystart(
        "OSA bug test",
        (("BAM", in_bam),),
        (("TSV", out_reads_info),),
        version=False,
    )
    verboselog = cli_utils.get_verboselog(fancylog, True)
    bf = pysam.AlignmentFile(in_bam, "rb")

    # Keeps track of the names of reads with OSAs.
    # We'll filter these reads out completely, so that they are not represented
    # in any of the alignments remaining in the BAM file.
    #
    # I was initially going to maintain a sorted array for this, so that we
    # could use binary search to quickly check if reads seen in the final
    # "pass" over the BAM file had OSAs, but it turns out that in Python using
    # sets is probably a better idea (or at the very least a simpler one):
    # https://stackoverflow.com/a/212971
    reads_with_osa = set()

    impacted_reads = set()

    # If literally nothing is aligned to this seq (before filtering), set its
    # entry in this to True -- this way, we can skip some work later
    seq2isempty = defaultdict(bool)

    with open(out_reads_info, "w") as f:
        f.write("Contig\tRead Name\tAlnCoords_0Indexed_Inclusive_Correct\n")

    for si, seq in enumerate(bf.references, 1):
        cli_utils.proglog(
            seq,
            si,
            bf.nreferences,
            verboselog,
            prefix="On ",
        )

        # Identify all linear alignments of each read to this sequence
        num_lin_alns = 0
        readname2Coords = defaultdict(list)
        readname2BadCoords = defaultdict(list)
        for num_lin_alns, linearaln in enumerate(bf.fetch(seq), 1):
            rn = linearaln.query_name

            # this uses the correct definition...!
            alncoords = bam_utils.get_coords(linearaln)
            readname2Coords[rn].append(alncoords)

            # (this doesn't)
            readname2BadCoords[rn].append((alncoords[0], alncoords[1] + 2))

        # How many (unique) reads are aligned total to this contig?
        n_reads_in_seq = len(readname2Coords)
        if n_reads_in_seq == 0:
            seq2isempty[seq] = True
            verboselog(
                f"Nothing is aligned to contig {seq}! Ignoring this contig.",
                prefix="",
            )
            continue

        # Sanity checking -- should never happen (TM) because we should have
        # already continued if n_reads_in_seq == 0
        if num_lin_alns == 0:
            raise WeirdError(
                "0 linear alns, but > 0 aligned reads? Something's wrong."
            )

        # Identify overlapping alignments from the same read
        n_reads_w_osa_in_seq = 0
        n_impacted_reads_in_seq = 0
        for rn in readname2Coords:
            alns = readname2Coords[rn]
            if len(alns) > 1:
                this_aln_really_has_osa = False
                # Okay, so this particular read has multiple supplementary
                # alignments to this sequence. Check if they overlap.
                # I have a feeling that examining all pairs of alignments is
                # an inefficient way of checking this, in theory -- however, in
                # practice most reads should probably have at most, like, 4
                # distinct linear alignments to a sequence, so this shouldn't
                # be a big bottleneck. Probably.
                for (a1, a2) in combinations(alns, 2):
                    # Efficiently test for overlap between two ranges:
                    # https://stackoverflow.com/a/3269471
                    if a1[0] <= a2[1] and a2[0] <= a1[1]:
                        # Okay, these two alignments of this read overlap.
                        reads_with_osa.add(rn)
                        n_reads_w_osa_in_seq += 1
                        this_aln_really_has_osa = True
                        break

                if not this_aln_really_has_osa:

                    # if we've made it here, then we know that this read does not
                    # have an OSA (using the correct definition). see if this read
                    # has an OSA using the incorrect, buggy definition (with
                    # off-by-2 end coordinates)
                    bad_alns = readname2BadCoords[rn]
                    for (a1, a2) in combinations(bad_alns, 2):
                        if a1[0] <= a2[1] and a2[0] <= a1[1]:
                            impacted_reads.add(rn)
                            with open(out_reads_info, "a") as f:
                                # output the *correct* alns
                                f.write(f"{seq}\t{rn}\t{alns}\n")
                            n_impacted_reads_in_seq += 1
                            break

        verboselog(
            f"There are {num_lin_alns:,} linear alignment(s) (from "
            f"{n_reads_in_seq:,} unique read(s)) to contig {seq}.",
            prefix="",
        )
        # We can compute this percentage without worrying about division by
        # zero because we've already ensured above that n_reads_in_seq != 0.
        rpct = 100 * (n_reads_w_osa_in_seq / n_reads_in_seq)
        verboselog(
            f"{n_reads_w_osa_in_seq:,} / {n_reads_in_seq:,} ({rpct:.2f}%) "
            "of these unique read(s) have OSAs.",
            prefix="",
        )
        bpct = 100 * (n_impacted_reads_in_seq / n_reads_in_seq)
        verboselog(
            f"{n_impacted_reads_in_seq:,} / {n_reads_in_seq:,} ({bpct:.2f}%) "
            "of these unique read(s) would have been incorrectly labelled as "
            "having OSAs.",
            prefix="",
        )

    bf.close()

# check on the sheepgut dataset
# check_osa_bug(
#     "/Poppy/mfedarko/sheepgut/sf-analyses/sheep/output/aln/sorted-unfiltered.bam",
#     "/Poppy/mfedarko/sheepgut/sf-analyses/sheep/output/osa-bug-read-info.tsv"
# )

# ... and on the chickengut dataset
check_osa_bug(
    "/Poppy/mfedarko/chicken-gut-meta/sf/aln-post-osa-bug-fix/sorted-unfiltered.bam",
    "/Poppy/mfedarko/sheepgut/sf-analyses/chicken/output/osa-bug-read-info.tsv"
)
