# SheepGut analysis scripts/notebooks

This repository contains various code files used in the analysis of metagenome-assembled genomes in a high-coverage sheep gut metagenome.

Note that many of these files will not work on other systems "out of the box"
(i.e. without some configuration). There are two main reasons for this.
First, we have not provided many of the input data files (which are well over
GitHub's 100 MB filesize limit as of writing) in this repository; and second,
many of these scripts assume that they are being run on our computing server,
and therefore use some things like hardcoded file paths that will break on
other systems.

In spite of this, we've tried to make things clear; if you have any questions
about getting things running or how things in this repository work, please feel
free to reach out.

## `main-workflow/`

This directory contains various scripts that work with the input data
(reads and an assembly graph) to produce statistics detailing the _mutation
spectrum_, as well as some other important values, throughout some genomes of
interest.

The basic "pipeline," which is performed by `RUN-ME.sh` in this directory,
involves the following general steps:

- Creating a FASTA file of all the edge sequences in a metagenome assembly graph

- Aligning the raw sequencing data (i.e. reads) back to the edge sequences in
  the metagenome assembly graph, creating a SAM file

- Filtering the resulting SAM file to remove secondary alignments (but not supplementary alignments!), and converting it to a BAM file

- Sorting and indexing the BAM file

- Filtering the BAM file to remove partially-mapped reads, and reads mapped completely or almost completely to other edge sequences we are not interested in

- Sorting and indexing the (filtered) BAM file

- Converting the final BAM file to a collection of simple
  [JSON](https://en.wikipedia.org/wiki/JSON) files representing the number
  of matches, mismatches, etc. at each position within the sequences of
  interest

The JSON files produced at the final step are the primary output from this
"main workflow."

## `inspect-seqs/`

This directory contains scripts that perform some tasks outside of the
"main workflow" involving analysis of the selected genomes/sequences
that don't really fit neatly anywhere else. Currently this includes scripts
for:

1. Running Prodigal to predict protein-coding genes within these sequences
2. Running CheckM

Note that, unlike the JSON files and other outputs of the "main workflow"
scripts above, these results are included in the repository since they're
both 1) relatively small files and 2) useful to have around in version control.

## `notebooks/`

This directory contains [Jupyter notebooks](https://en.wikipedia.org/wiki/Project_Jupyter#Jupyter_Notebook) that create various plots, files, etc. of the data.

Figures are output to a directory named `notebooks/figs/`.
Mutation profile tables are output to a directory named
`notebooks/mutation-profiles/`.

Many of these notebooks rely entirely or almost entirely on the JSON files
created by the data processing scripts above. The initial structure of this
repository was as two separate codebases, where I would create the JSONs on a
computing server and then download them to my laptop for further analysis;
however, eventually the JSONs became too large to load in memory on my poor
laptop, so I elected to run the notebooks on the server as well.

For reference, a script that runs these notebooks is located at `notebooks/RUN-NOTEBOOKS.py`.

## `find-edges-in-other-graph/`

This directory contains some code used in the transition from working with
a smaller assembly graph to working with a larger assembly graph (constructed
using the original data plus ~5x as much data).

One of the first tasks in this transition was taking sequences of interest
from the smaller graph and trying to find them in the larger graph -- code
for this is included here.
