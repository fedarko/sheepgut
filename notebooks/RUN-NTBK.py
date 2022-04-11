#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)


import sys
import time
import subprocess

ntbk_to_run = sys.argv[1]

print(f"Running notebook {ntbk_to_run}...")

t0 = time.time()

subprocess.run([
    "jupyter",
    "nbconvert",
    "--execute",
    "--ExecutePreprocessor.timeout=None",
    "--ExecutePreprocessor.kernel_name=python3",
    "--to",
    "notebook",
    "--inplace",
    ntbk_to_run
], check=True)

t1 = time.time()
print("Done. Took {t1 - t0:,.2f} seconds.")
