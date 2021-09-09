#! /usr/bin/env python3

import time
import glob
import skbio

OUTFILE = "../notebooks/misc-text/read-stats.tex"

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

        read_kwargs = {}
        # Very primitive "sniffer" for the filetype. This is a very lazy way
        # to do this, since misnamed files (e.g. a FASTQ that's named
        # evil-file.fasta.gz) will cause problems. So this should only be used
        # with correctly named files! ... However, skbio.io.read()'s default
        # use of verify=True should mean it'll throw an error if weird stuff
        # is detected, so we're decently safe.
        irf_lc = irf.lower()

        if irf_lc.endswith("fastq.gz"):
            ft = "fastq"
            # skbio.io.read() needs us to pass either the "variant" or
            # "phred_offset" argument when reading FASTQ files. We don't care
            # about quality scores here (at least in the case of counting read
            # lengths), so we just set "variant" to "sanger" since skbio
            # doesn't (to my knowledge) have a specific variant name set up for
            # PacBio HiFi / CCS.
            read_kwargs["variant"] = "sanger"

        elif irf_lc.endswith("fasta.gz"):
            ft = "fasta"

        else:
            raise ValueError("Can't tell what filetype this file is?")

        file_num_reads = 0
        file_read_length = 0
        # Use of skbio.io.read() to get a generator derived from
        # github.com/biocore/scikit-bio/issues/1418#issuecomment-241498418
        # Unpacking of keyword arguments here derived from
        # https://stackoverflow.com/a/50258379
        read_gen = skbio.io.read(
            irf, constructor=skbio.DNA, format=ft, **read_kwargs
        )
        for read_num, read in enumerate(read_gen, 1):
            if read_num == 1 or read_num % 1000 == 0:
                print(f"On read {read_num:,} in file #{num_files_seen}...")
            file_num_reads += 1
            file_read_length += len(read)

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

print(f"Writing out results to {OUTFILE}...", end=" ", flush=True)
with open(OUTFILE, "w") as of:
    of.write(
        f"In total, these HiFi data include {total_num_reads:,} reads\n"
        f"with total length {total_read_length:,} bp\n"
        f"and average length {avg_read_length:,.2f} bp."
    )
print("Done!")
