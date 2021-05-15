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
# -Notebooks that rely on Graphviz' sfdp utility and the PRISM overlap
#  removal method -- which in turn relies on the "gts" library -- which I
#  have not been able to get installed on our cluster yet, so I've just been
#  running this notebook on my laptop for now (linked mutation call&plot)
# -Notebooks that are only here so I can test stuff locally for convenience's
#  sake rather than pushing things up to the cluster (useful when messing with
#  matplotlib) (the "local version" of the graph notebook)
EXCLUDED_NOTEBOOKS = [
    "Header.ipynb",
    "LoadMutationJSONData.ipynb",
    "GeneUtils.ipynb",
    "Codon_AminoAcid_MutationMatrices.ipynb",
    "LinkedMutations-ClassifyReads.ipynb",
    "LinkedMutations-CallAndPlot.ipynb",
    "GraphNtbk_LocalVsn_ForTestingCovLenSummaryPlot.ipynb",
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
