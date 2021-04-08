# SheepGut Analysis Scripts/Notebooks

This repository contains various code files used in the analysis of metagenome-assembled genomes in a high-coverage sheep gut metagenome.

Note that many of these files will not work on other systems since we have not
provided the raw data files (which are over 100 GB) in this repository.
However, we've tried to document things clearly; if you have any questions,
please feel free to reach out.

## `data-processing-scripts/`

This directory contains various scripts that work with the raw data.
The basic workflow used involves the following steps:

- Creating a FASTA file of all the edge sequences in a metagenome assembly graph

- Aligning the raw sequencing data (i.e. reads) back to the edge sequences in
  the metagenome assembly graph

- Filtering the resulting alignment to remove:
  - Secondary alignments
  - Reads mapped to the edge sequences we are interested in that are only
    partially mapped to a particular edge sequence of interest (or to another
    edge in the same component of the assembly graph, if a sequence we are
    interested in is not in its own isolated component of the assembly graph)

- Sorting and indexing the alignment

- Converting the alignment to a collection of simple
  [JSON](https://en.wikipedia.org/wiki/JSON) files representing the number
  of matches, mismatches, etc. at each position within the sequences of
  interest

The JSON files are, as of writing, the primary output from these scripts.
However, there are other scripts included within this directory that generate
some other data (e.g. running Prodigal; running CheckM; analyzing certain
regions of these sequences as needed for the paper).

## `analysis-notebooks/`

This directory contains [Jupyter notebooks](https://en.wikipedia.org/wiki/Project_Jupyter#Jupyter_Notebook) that create various plots, files, etc. of the data.

Many of these notebooks rely entirely or almost entirely on the JSON files
created by the data processing scripts above. The initial structure of this
repository was as two separate codebases, where I would create the JSONs on a
computing server and then download them to my laptop for further analysis;
however, eventually the JSONs became too large to load in memory on my poor
laptop, so I elected to run the notebooks on the server as well.

## `find-edges-in-other-graph/`

This directory contains some code used in the transition from working with
a smaller assembly graph to working with a larger assembly graph (constructed
using the original data plus ~5x as much data).

One of the first tasks in this transition was taking sequences of interest
from the smaller graph and trying to find them in the larger graph -- code
for this is included here.
