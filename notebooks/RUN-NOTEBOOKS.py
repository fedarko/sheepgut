#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)
# NOTE: currently this just runs the notebooks in sorted filename order, which
# is good enough. Ideally we'd set this up using Snakemake / etc., since some
# notebooks rely on others having already been run (e.g. matrices, phasing).


import os
import subprocess

EXCLUDED_NOTEBOOKS = [
    "Header.ipynb",
    "GeneUtils.ipynb",
]

cwd_files = os.listdir()
# Make order deterministic by sorting
for f in sorted(cwd_files):
    if f.endswith(".ipynb"):
        if f.startswith("Phasing-"):
            print("Running notebook {}".format(f))
            subprocess.run([
                "jupyter",
                "nbconvert",
                "--execute",
                "--ExecutePreprocessor.timeout=None",
                "--ExecutePreprocessor.kernel_name=python3",
                "--to",
                "notebook",
                "--inplace",
                f
            ], check=True)
