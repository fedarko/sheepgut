#! /usr/bin/env python3
# (Sorry, this is copied from the ../sheep/ directory... ideally this would be
# abstracted to a single script but time is tight and i don't wanna abstract
# this and then find out that it broke the original use case in some small way)

import pysam
import pysamstats


BAM_FP = "/Poppy/mfedarko/chicken-gut-meta/sf/aln/final.bam"
FASTA_FP = "/Poppy/mfedarko/chicken-gut-meta/sf/contigs.fasta"

OUT_READ_INFO_FP = "output/deletion-data/read-info.tsv"
OUT_CONTIG_INFO_FP = "output/deletion-data/contig-info.tsv"

# ... everything below this line is the same as the version i copied over from
# ../sheep/, as of writing at least

# Initialize read TSV file
with open(OUT_READ_INFO_FP, "w") as f:
    f.write("ReadName\tContig\tRefStart\tRefEnd\tNumDels\n")
    
# Initialize contig TSV file

# count all positions in a contig with >= 5 deletions, >= 10 deletions, ...
# THIS MUST BE IN SORTED ASCENDING ORDER or else it'll break stuff below
del_thresholds = [5, 10, 25, 50, 100, 200, 500, 1000]
assert del_thresholds == sorted(del_thresholds) and len(set(del_thresholds)) == len(del_thresholds)

with open(OUT_CONTIG_INFO_FP, "w") as f:
    header = "Contig"
    for d in del_thresholds:
        header += f"\t{d}DelPos"
    header += "\n"
    f.write(header)

# OK now let's rip through this BAM file faster than I just ripped into one of
# those bags of microwave popcorn (it's "microwave" popcorn, right? I feel like
# if I say "microwaveable" Merriam and Webster will pop out from behind the
# door and punch me in the face)

bf = pysam.AlignmentFile(BAM_FP)

for ci, contig in enumerate(bf.references, 1):
    print(f"On contig {contig} ({ci:,} / {bf.nreferences:,})...")
    
    read_tsv_txt = ""
    for ai, aln in enumerate(bf.fetch(contig), 1):
        num_dels = 0
        for op_len_pair in aln.cigartuples:
            if op_len_pair[0] == pysam.CDEL:
                # this is a deletion operation: see
                # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
                # (We could also just say 2, since pysam.CDEL == 2, but this way is clearer and
                # probs less prone to breaking horrifically lol)
                num_dels += op_len_pair[1]
        
        # this uses aln.reference_length, which is just based on the aligned parts of the read in this
        # linear alignment -- *not* the full read length. i guess we could also do this the other way tho
        read_tsv_txt += f"{aln.query_name}\t{contig}\t{aln.reference_start}\t{aln.reference_end}\t{num_dels}\n"
        
    with open(OUT_READ_INFO_FP, "a") as f:
        f.write(read_tsv_txt)
        
    # max depth here matches sf config as of writing:
    # https://github.com/fedarko/strainFlye/blob/main/strainflye/config.py
    del_cts = [0] * len(del_thresholds)
    for pos, rec in enumerate(
        pysamstats.stat_variation(bf, chrom=contig, fafile=FASTA_FP, pad=True, max_depth=100000000)
    ):
        pos_del_ct = rec["deletions"]
        for di, d in enumerate(del_thresholds):
            if pos_del_ct >= d:
                del_cts[di] += 1
            else:
                # can break early, since we sorted del_thresholds in ascending order
                break

    with open(OUT_CONTIG_INFO_FP, "a") as f:
        row_txt = contig
        for d in del_cts:
            row_txt += f"\t{d}"
        row_txt += "\n"
        f.write(row_txt)
