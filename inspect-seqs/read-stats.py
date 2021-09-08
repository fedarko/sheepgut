#! /usr/bin/env python3

import time
import glob
import skbio

t0 = time.time()
print("Starting processing files...")
with open("../config/input-reads", "r") as input_reads_filepaths_file:
    ir_paths = input_reads_filepaths_file.readline().strip().split()

num_files_seen = 0
total_num_reads = 0
total_read_length = 0

# Go through each general path -- this could be a single file, or a reference
# to multiple files (e.g. "myfiles/*.fastq.gz", which could correspond to
# arbitrarily many files)
for irp in ir_paths:
    # Use glob.glob() to convert these to lists of real file paths.
    # This is of course vulnerable to a race condition if some antagonist
    # deletes a file in between calling glob.glob() here and when we attempt to
    # parse it below. But that should trigger an error, and ... look, this is a
    # script for counting read lengths for a single paper. It should be fine.
    for irf in glob.glob(irp):
        tf0 = time.time()
        num_files_seen += 1
        print(f"On file {irf} (file #{num_files_seen:,})...")

        file_num_reads = 0
        file_read_length = 0
        # github.com/biocore/scikit-bio/issues/1418#issuecomment-241498418
        seq_gen = skbio.io.read(irf, constructor=skbio.DNA)
        for seq in seq_gen:
            file_num_reads += 1
            file_read_length += len(seq)

        tf1 = time.time()
        print(f"Done. Processing this file took {tf1 - tf0:,.2f} sec.")
        print(
            f"This file had {file_num_reads:,} reads of total length "
            f"{file_read_length:,}."
        )
        total_num_reads += file_num_reads
        total_read_length += file_read_length

avg_read_length = total_read_length / total_num_reads

t1 = time.time()
print(f"Done: processed {num_files_seen:,} files.")
print(f"Total time taken: {t1 - t0:,.2f} sec.")

print("=" * 79)
print(f"Total number of reads: {total_num_reads:,}")
print(f"Total read length: {total_read_length:,}")
print(f"Average read length: {avg_read_length:,}")
