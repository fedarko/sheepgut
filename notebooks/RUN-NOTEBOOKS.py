#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)

import os
import subprocess

# Notebooks we don't (explicitly) run here. Many of these are notebooks
# that exist only to be run by other notebooks (e.g. header, gene utils,
# load mutation JSON data); this could also be used for unfinished notebooks,
# stuff that I don't need to frequently rerun, etc.
EXCLUDED_NOTEBOOKS = [
    "Header.ipynb",
    "GraphCoverageAndConnectivityAnalysis.ipynb",
    "LoadMutationJSONData.ipynb",
    "GeneUtils.ipynb",
    "Codon_AminoAcid_MutationMatrices.ipynb",
]

cwd_files = os.listdir()
for f in cwd_files:
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
            ])
