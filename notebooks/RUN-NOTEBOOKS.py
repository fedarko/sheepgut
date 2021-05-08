#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)

import os
import subprocess

# Notebooks we don't (explicitly) run here. This includes:
# -Notebooks that exist only to be run by other notebooks
#  (e.g. header, gene utils, load mutation JSON data)
# -Notebooks that take a super long time to rerun, and should therefore be
#  rerun only when definitely needed
#  (e.g. codon/aa mutation matrices, linked mutation read classification)
EXCLUDED_NOTEBOOKS = [
    "Header.ipynb",
    "LoadMutationJSONData.ipynb",
    "GeneUtils.ipynb",
    "Codon_AminoAcid_MutationMatrices.ipynb",
    "LinkedMutations-ClassifyReads.ipynb",
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
