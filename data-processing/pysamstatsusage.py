import pysam
import pysamstats
bf =  pysam.AlignmentFile("aln-sorted.bam")
bf
for rec in pysamstats.stat_coverage(bf, chrom="edge_6104", start=200, end=210):
    print(rec)
for rec in pysamstats.stat_coverage(bf, chrom="edge_6104", start=200, end=210, truncate=True):
    print(rec)
for rec in pysamstats.stat_coverage(bf, chrom="edge_6104", start=200, end=210, truncate=True, no_del=True):
    print(rec)
for rec in pysamstats.stat_coverage(bf, chrom="edge_6104", start=200, end=210, truncate=True, no_deltaco=True):
    print(rec)
help(pysamstats.stat_coverage)
help(pysamstats)
for rec in pysamstats.stat_variation(bf, chrom="edge_6104", start=200, end=210, truncate=True):
    print(rec)
for rec in pysamstats.stat_variation(bf, chrom="edge_6104", fafile="newseqs.fasta", start=200, end=210, truncate=True):
    print(rec)
%history 1-13 -f pysamstatsusage.py
