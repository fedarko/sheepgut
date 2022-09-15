# SheepGut analysis scripts/notebooks

[![SheepGut Analysis Code CI](https://github.com/fedarko/sheepgut/actions/workflows/main.yml/badge.svg)](https://github.com/fedarko/sheepgut/actions/workflows/main.yml)

This repository contains various code files used in the analysis of metagenome-assembled genomes in a high-coverage sheep gut metagenome.

Broadly, the "inputs" to these analyses are:

1. An assembly graph (in GFA 1 format)
2. The reads used to generate the assembly graph (in FASTQ and/or FASTA format)

**Please note** that many of the analyses shown here (alignment, mutation
calling, FDR estimation, hotspot / coldspot finding, smoothed read creation and
assembly) have been ported to the
[strainFlye](https://github.com/fedarko/strainFlye) repository.
strainFlye's code should be easier to use and more well-tested than the code
here -- so I recommend using strainFlye instead of this messy code if possible.
(These instructions are primarily here for the sake of reproducibility, and
because as of writing the code for
link graphs + mutation matrices + growth dynamics are not ported over to
strainFlye yet.)

## Reproducing these analyses on your own system

If you have any questions about getting things running or how things in this
repository work, please feel free to reach out. This code was primarily run on
an Ubuntu 16.04 system.

### Installation

Most of the relevant "dependency" software
for these analyses we used (e.g. minimap2,
samtools, pysam, ...) was installed using [conda](https://conda.io/) (and
[pip](https://pip.pypa.io/) inside of conda, in some cases). This section
assumes you have conda installed.

First, let's download this repository's code:

```bash
git clone https://github.com/fedarko/sheepgut.git
cd sheepgut
```

Now, we need to install dependencies. We use
[mamba](https://mamba.readthedocs.io/en/latest/index.html) to do this because
it's fast, but using the default conda installer should also work fine.

```bash
# Install mamba: https://mamba.readthedocs.io/en/latest/installation.html
conda install mamba -n base -c conda-forge

# Now install our packages using mamba
mamba env create -f environment.yml
```

This creates a conda environment named `sheepgut` that you can use for the bulk
of these analyses. (Exceptions -- that is, analyses that require additional
software to be installed -- are detailed below.)

### Downloading input data and updating the `config` folder

The main two inputs to these analyses are reads and an assembly graph
produced from these reads.
You can download a copy of the assembly graph from
[Zenodo](https://zenodo.org/record/6545142). Once you've downloaded it, please
update the corresponding path in `config/input-graph` (a few of the scripts
will make use of this path).

Similarly, you should then download reads and update the path in
`config/input-reads` accordingly. (Since the reads will probably span multiple
files, you can separate filenames by spaces -- this is done in
`config/input-reads`, as of writing.)

### Alignment and alignment filtering

Run the following commands in the terminal:

```bash
# Run "main workflow" (do alignment, etc.)
cd main-workflow/
./RUN-ME.sh
```

This will create a folder in `main-workflow/` named `output/`; this folder will
contain a BAM file corresponding to the filtered alignment
(`main-workflow/output/fully-filtered-and-sorted-aln.bam`), as well as a
"pickle" file containing pileup data for the three selected MAGs
(`main-workflow/output/seq2pos2pileup.pickle`).

(The "pickle" file is an inefficient way of making pileup access easy; the
actual strainFlye codebase doesn't need to use this, since it just uses pysam
directly.)

### Run various read / MAG-level analyses

These steps are optional (we include CheckM, Prodigal, barrnap, LoFreq, and
read stats outputs in this repository). We assume that you start each command
from the root of the repository.

#### CheckM (compute completeness and contamination statistics for the three MAGs)

```bash
cd inspect-seqs/
./checkm.sh
```

#### Prodigal (predict protein-coding genes in the three MAGs)

```bash
cd inspect-seqs/
./prodigal.py
```

#### Barrnap (predict rRNA genes in BACT1)

```bash
cd inspect-seqs/
./barrnap.sh
```

#### LoFreq

Note that you will need to install [LoFreq](http://csb5.github.io/lofreq/) to
run this script, and that this is not included in the `environment.yml` file
shown above. When I needed to run LoFreq, I wound up creating a new conda
environment just for installing / running LoFreq.

(Also, this assumes that you have already run `main-workflow/RUN-ME.sh` and
created an alignment.)

```bash
cd main-workflow/

# Create a BAM file for just the three selected MAGs (saves time)
./filter-bam-to-selected-mags.sh

# Run LoFreq using it
./run-lofreq.py
```

#### Read statistics

```bash
cd inspect-seqs/
./read-stats.py
```

### Run analysis notebooks

We have now set the stage for running the analysis notebooks in this repository.
These notebooks contain most of the "higher-level" analyses shown in our paper.

You can run notebooks in the `notebooks/` folder by starting a Jupyter notebook server
(see e.g. [this documentation](https://docs.jupyter.org/en/latest/running.html))
from within this directory, then opening up any of the notebooks. You can also
run individual notebooks from the command line directly
using the `RUN-NTBK.py` Python script located within this repository (this may
be useful for some of the notebooks which take a while to run).

#### Notebook outputs

Most figures are output to a directory named `notebooks/figs/`, although there are a
few other types of outputs created. For example, some notebooks output things
to the `notebooks/misc-text/` directory, which contains `.tex` files that are
included in our paper.

#### Notebooks that depend on other notebooks already having been run

Please note that some notebooks depend on the outputs of other notebooks having
already been created. I have tried to track these "dependencies" here. (I have
tried to include the "intermediate" outputs here in the repository, when
possible; however, the `notebooks/phasing-data/` outputs are not included
because they take up a lot of space.)

- `DemonstratingTargetDecoyApproach.ipynb`
  - Depends on the output of
    `SynAndNonsenseMutationRateBarplots-PositionBased.ipynb` (the
    `notebooks/misc-output/pos_p2seq2*.pickle` files).

- `Matrices-02-Viz.ipynb` 
  - Depends on `Matrices-01-Compute.ipynb` (the `notebooks/matrix-jsons/*`
    files).

- `Phasing-02-VizGraph.ipynb`
  - Depends on `Phasing-01-MakeGraph.ipynb` (the
    `notebooks/phasing-data/*.pickle` files).

- `Phasing-03-MiscStats.ipynb`
  - Depends on `Phasing-01-MakeGraph.ipynb` (the
    `notebooks/phasing-data/*.pickle` files).

- `Phasing-LJA-SmoothedVirtualReadCoverage.ipynb` (unused in the paper, as of
  writing)
  - Depends on `Phasing-LJA.ipynb` (the
    `phasing-data/smoothed-reads/edge_1671_smoothed_reads_delignore_vrlow.fasta`
    file).

#### Notebooks that need extra software to be installed

- `notebooks/Phasing-02-VizGraph.ipynb`
  - Depends on [Graphviz](https://www.graphviz.org/) -- in particular, the
    `sfdp` layout algorithm -- being installed. (Note that the part of this
    notebook that actually runs `sfdp` is commented out at the moment; so, for
    now, this notebook will just generate `.gv` files without actually
    visualizing them.)

- `notebooks/Phasing-LJA.ipynb` and `notebooks/Phasing-LJA-CAMP-Gene.ipynb`
  - Depend on [LJA and jumboDBG](https://github.com/AntonBankevich/LJA/)
    being installed.
