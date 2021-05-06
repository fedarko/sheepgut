# SheepGut analysis scripts/notebooks

This repository contains various code files used in the analysis of metagenome-assembled genomes in a high-coverage sheep gut metagenome.

Broadly, the "inputs" to these analyses are:

1. An assembly graph (in GFA format)
2. The reads used to generate the assembly graph (in FASTQ and/or FASTA format)

## Reproducing these analyses on your own system

If you have any questions about getting things running or how things in this
repository work, please feel free to reach out. This guide should work on most
Linux / macOS systems; it's been tested on Ubuntu 16.04.

### Setting up your environment
Most of the relevant software for these analyses we used (e.g. minimap2,
samtools, pysam, ...) was installed using [conda](https://conda.io/) (and
[pip](https://pip.pypa.io/) inside of conda, in some cases).

The `conda-environment.yml` file in this repository details the conda
environment that the analyses were ran within.
(This file was created using [`conda env export`](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment).)

You should be able to replicate this environment from the YML file directly
(see [the conda docs on how to do this](https://docs.conda.io/projects/conda/en/master/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)).

**NOTE: The conda environment in the YML file is pretty bloated -- there are a
few things in there that are not needed for these analyses. This envirnoment
should be replaced with a smaller, minimal one.**

You could also create a fresh Python 3.7 environment and install specific
needed packages individually. However, that might be tedious and prone to
errors.

### Basic walkthrough

1. Download assembly graph; update the path in `config/input-graph` accordingly
2. Download reads; update the path in `config/input-reads` accordingly
3. Run the following commands in bash:

```bash
# Run "main workflow" (do alignment, etc.)
cd main-workflow/
./RUN-ME.sh

# Run various sequence-level analyses
cd ../inspect-seqs/
./checkm.sh
./prodigal.py

# Create figures for the paper
cd ../notebooks
./RUN-NOTEBOOKS.py
```

### Input file locations

Note that we have not provided the input data files (or many of the
"intermediate" data files) in this repository due to GitHub's filesize limits.
You'll need to download them yourself and update some paths accordingly (see
above).

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
