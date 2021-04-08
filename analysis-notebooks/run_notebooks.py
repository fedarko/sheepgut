#! /usr/bin/env python3
# From https://github.com/biocore/qurro/blob/master/Makefile
# (which in turn was from https://nbconvert.readthedocs.io/en/latest/usage.html#notebook-and-preprocessors)

import os
import subprocess

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
