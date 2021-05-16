#! /usr/bin/env python3

import glob
import subprocess

with open("../config/input-reads", "r") as rfc:
    read_files_unexpanded = next(rfc).strip().split()

total_num_reads = 0
total_read_length = 0
# For each unexpanded read file path...
for rfu in read_files_unexpanded:
    # Expand this path, and iterate through all of the corresponding files!
    for rf in glob.glob(rfu):
        print(f"On file {rf}...")
        # Now we're looking at a single one of the read files.
        num_lines_s = subprocess.run(
            f"zcat {rf} | wc -l",
            capture_output=True,
            shell=True
        )
        num_lines = int(num_reads_s.strip())
        if ".fastq" in rf:
            if num_lines % 4 != 0:
                raise ValueError(
                    f"File {rf} has {num_lines} lines; not divisible by 4???"
                )
            num_reads = num_lines / 4
        else:
            num_reads = num_lines / 2
