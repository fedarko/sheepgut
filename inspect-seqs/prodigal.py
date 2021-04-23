#! /usr/bin/env python3

# running with the default normal (in this version, "single") mode since these
# are long enough sequences. We're splitting them up into multiple runs per the
# Prodigal wiki's advice --
# https://github.com/hyattpd/prodigal/wiki/advice-by-input-type.

import os
import subprocess

# os.scandir() is a slightly fancier way of iterating over files in a directory
# than os.listdir() since it makes it easier to e.g. distinguish files from
# directories.
for de in os.scandir("../seqs"):
    if de.is_file():
        fn = de.name
        if fn.lower().endswith(".fasta"):
            in_file_path = os.path.join("..", "seqs", fn)
            # We use the [:-5] to trim off the ".fasta" to make using
            # different file suffixes easier
            out_file_path_base = os.path.join(
                "..", "seqs", "genes", fn[:-6]
            )
            print("Running Prodigal on sequence {}".format(fn))
            subprocess.run([
                "prodigal",
                "-i",
                in_file_path,
                "-o",
                out_file_path_base + ".sco",
                "-a",
                out_file_path_base + "_aa.fasta",
                "-f",
                "sco"
            ])
        else:
            print(
                (
                    "WARNING: Found file {} in seqs/ without .fasta extension?"
                    "Not running Prodigal on it, at least."
                ).format(fn)
            )
