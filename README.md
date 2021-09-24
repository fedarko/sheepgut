# SheepGut analysis scripts/notebooks

[![SheepGut Analysis Code CI](https://github.com/fedarko/sheepgut/actions/workflows/main.yml/badge.svg)](https://github.com/fedarko/sheepgut/actions/workflows/main.yml)

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
[pip](https://pip.pypa.io/) inside of conda, in some cases). This section
assumes you have at least conda installed.

It will probably be fastest and easiest to install the dependencies for this
project using conda via the following command.
This command will not install the exact same versions of these
packages as I've used in development, but the versions installed should be
close enough:

```bash
conda env create -f environment.yml
```

If you'd prefer to try to use the _exact_ same package versions I had
installed, this environment contains the exact versions of the necessary
dependencies. Please note that this "exact" conda environment file is pretty
bloated; there are many things in there that are not needed for these analyses.

```bash
conda env create -f exact-environment.yml
```

I am aware that it would probably be best to set up a Docker container or
something -- getting to that is on my radar, but I am not sure I'll have time.

### Basic walkthrough

This will walk you through reproducing the results after setting up your
environment.

1. Download assembly graph; update the path in `config/input-graph` accordingly
2. Download reads; update the path in `config/input-reads` accordingly
3. Run the following commands in bash:

```bash
# Run "main workflow" (do alignment, etc.)
cd main-workflow/
./RUN-ME.sh

# Run various read / sequence-level analyses
cd ../inspect-seqs/
./read-stats.py # compute read length statistics
./checkm.sh     # Run CheckM on certain MAGs
./prodigal.py   # Run Prodigal on certain MAGs

# Create figures for the paper
# (Note that this won't run all of the notebooks, since some of them
# take a relatively long time to run; see the script's comments for details)
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
"main workflow" involving analysis of the sequencing data
that don't really fit neatly anywhere else. Currently this includes scripts
for:

1. Computing read length statistics
2. Running CheckM
3. Running Prodigal to predict protein-coding genes within these sequences

Note that, unlike the JSON files and other outputs of the "main workflow"
scripts above, the results from running CheckM and Prodigal are
included in the repository since they're
both 1) relatively small files and 2) useful to have around in version control.
(The read length stats script doesn't create any files; it just outputs
information to stdout.)

## `notebooks/`

This directory contains [Jupyter Notebooks](https://en.wikipedia.org/wiki/Project_Jupyter#Jupyter_Notebook) that create various plots, files, etc. of the data.

Figures are output to a directory named `notebooks/figs/`.
Mutation profile tables are output to a directory named
`notebooks/mutation-profiles/`.

Many of these notebooks rely entirely or almost entirely on the JSON files
created by the data processing scripts above. The initial structure of this
repository was as two separate codebases, where I would create the JSONs on a
computing server and then download them to my laptop for further analysis;
however, eventually the JSONs became too large to load in memory on my poor
laptop, so I elected to run the notebooks on the server as well.

For reference, a script that runs some of these
notebooks is located at `notebooks/RUN-NOTEBOOKS.py`.

## `find-edges-in-other-graph/`

This directory contains some code used in the transition from working with
a smaller assembly graph to working with a larger assembly graph (constructed
using the original data plus ~5x as much data).

One of the first tasks in this transition was taking sequences of interest
from the smaller graph and trying to find them in the larger graph -- code
for this is included here.
