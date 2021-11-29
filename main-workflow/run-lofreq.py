#! /usr/bin/env python3
# Runs LoFreq on just the three selected MAGs, while timing it.
#
# NOTE 1: I wound up installing LoFreq in a separate conda environment than
# the big one I use for most of the other analyses here, so this is a reminder
# to future-me to activate that environment before running this script.
# NOTE 2: LoFreq will be "cowardly" (its words, not mine :) and refuse to
# overwite the output VCF file, so you've gotta manually remove that before
# re-running this script.
import time
import subprocess

# Using the timeit() module would be ideal, but running LoFreq just once
# already takes quite a lot of time.
t0 = time.time()
print(f"Start time: {t0:,} sec.")
print("Running LoFreq on the three selected MAGs...")
subprocess.run(
    [
        "lofreq",
        "call",
        "-f",
        "output/selected-mags.fasta",
        "-o",
        "output/lofreq.vcf",
        "output/selected-mags.bam"
    ],
    check=True
)
t1 = time.time()
print(f"End time: {t1:,} sec.")
tt = t1 - t0
print(f"Time taken: {tt:,.2f} seconds.")
with open("../notebooks/misc-text/lofreq-runtime-3-mags.tex", "w") as tf:
    tf.write(f"{tt:,.2f} seconds\\endinput")
