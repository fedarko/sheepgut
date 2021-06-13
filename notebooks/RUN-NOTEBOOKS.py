#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)
# NOTE: currently this just runs the notebooks in sorted filename order, which
# is good enough. Ideally we'd set this up using Snakemake / etc., since some
# notebooks rely on others having already been run
# (so far the only instances of this I can remember are
# LinkedMutations-ClassifyReads before LinkedMutations-CallAndPlot, and
# MutationMatrices before PlotMutationMatrices).


import os
import subprocess

# Notebooks we don't (explicitly) run here. This includes:
# -Notebooks that exist only to be run by other notebooks
#  (e.g. header, gene utils, load mutation JSON data)
# -Notebooks that rely on Graphviz' sfdp utility and the PRISM overlap
#  removal method -- which in turn relies on the "gts" library -- which I
#  have not been able to get installed on our cluster yet, so I've just been
#  running this notebook on my laptop for now (linked mutation call&plot)
EXCLUDED_NOTEBOOKS = [
    "Header.ipynb",
    "LoadMutationJSONData.ipynb",
    "GeneUtils.ipynb",
    "LinkedMutations-CallAndPlot.ipynb",
]

cwd_files = os.listdir()
# Make order deterministic by sorting
for f in sorted(cwd_files):
    if f.endswith(".ipynb"):
        if f not in EXCLUDED_NOTEBOOKS:
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
