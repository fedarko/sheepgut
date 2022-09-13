#! /usr/bin/env python3

# This assumes that you're running this in the same directory as the output of
# "strainFlye call". That, or you can modify the "diversity-indices.tsv" string
# to point to the filepath to a diversity indices TSV file.
import pandas as pd
from statistics import mean
df = pd.read_csv("diversity-indices.tsv", sep="\t", index_col=0)
long = df[df["Length"] >= 1000000].sort_values(["AverageCoverage"], ascending=False)
# how many "long" contigs?
len(long)
# show how many long contigs have lengths above each threshold
long[long["AverageCoverage"] >= 200]
long[long["AverageCoverage"] >= 150]
long[long["AverageCoverage"] >= 100]
# mean coverage of long contigs
mean(long["AverageCoverage"])
