# strainFlye paper analysis scripts/notebooks

[![Analysis Code CI](https://github.com/fedarko/sheepgut/actions/workflows/main.yml/badge.svg)](https://github.com/fedarko/sheepgut/actions/workflows/main.yml)

This repository contains various code and data files used in the strainFlye
paper.

## Datasets

Most of these analyses are focused on the analysis of a sheep gut
metagenome dataset
([Kolmogorov _et al._, 2020](https://www.nature.com/articles/s41592-020-00971-x),
[Bickhart/Kolmogorov _et al._, 2022](https://www.nature.com/articles/s41587-021-01130-z));
however, some analyses (in particular, the files in `sf-analyses/chicken/`)
focus on a later analysis we performed of a chicken gut metagenome dataset
([Feng _et al._, 2022](https://www.nature.com/articles/s41592-022-01478-3)).

## The strainFlye pipeline

**Please note** that many of the analyses shown here (alignment, mutation
calling, FDR estimation, hotspot / coldspot finding, smoothed read creation and
assembly) have been ported to the
[strainFlye](https://github.com/fedarko/strainFlye) repository.

strainFlye's code should be easier to use and more well-tested than the code
here -- so, if you would like to replicate these sorts of analyses on another
dataset, I recommend using strainFlye instead of this messy code if possible.

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

Please note that, although the `config` folder is referenced by many of the
analyses in this repository, there are still many analyses here that directly
point to these or other SheepGut-related files. To reiterate the point above, I
strongly recommend using the
[strainFlye pipeline code](https://github.com/fedarko/strainFlye) for arbitrary
datasets.

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

The "pickle" file is an inefficient way of making pileup access easy for a
small set of MAGs; the actual strainFlye codebase doesn't need to use this,
since it just uses pysam directly.

### Run various read / MAG-level analyses

These steps are optional (we include CheckM, Prodigal, barrnap, read stats, and
LoFreq outputs in this repository). For each of these sets of commands, we
assume that you are located in the root of the repository.

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

#### Read statistics (output information about read lengths, etc.)

```bash
cd inspect-seqs/
./read-stats.py
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

# Output info on how many positions LoFreq called multiple variants at
# in the three selected MAGs
cd ../misc-scripts
./find_multivariant_pos_lofreq.py
```

### Run analysis notebooks

We have now set the stage for running the analysis notebooks in this repository.
These notebooks contain most of the "higher-level" analyses shown in our paper,
and create most of the figures shown there.

You can run notebooks in the `notebooks/` folder by starting a Jupyter notebook server
(see e.g. [this documentation](https://docs.jupyter.org/en/latest/running.html))
from within this directory, then opening up any of the notebooks. You can also
run individual notebooks from the command line directly
using the `RUN-NTBK.py` Python script located within this repository (this may
be useful for some of the notebooks which take a while to run).

#### Troubleshooting an "internal server error" when running Jupyter notebooks
When trying to log in to a Jupyter notebook session, you get an error that
shows up in the browser as `500: Internal Server Error`.
If the Jupyter server program produces an error that
looks like `ImportError: libffi.so.7: cannot open shared object file: No such file or
directory`, then this is a known issue with using conda and Jupyter -- see
[this GitHub issues
comment](https://github.com/conda/conda/issues/9038#issuecomment-627698375)
for details.

I ran into this error a few times. The solution resembles that described in the
GitHub issues comment above:

1. `cd` into the `lib/` folder for your conda environment: for an environment
   named `sheepgut`, for example, this is (at least on my system's conda
   installation) `cd ~/miniconda3/envs/sheepgut/lib/`

2. Figure out what versions of `libffi.so*` you have in this folder. Something
   like `ls -ahlt libffi.so*` should work -- at least for my installation, I had
   multiple versions of `libffi.so*` in this folder. However, all but one of
   these were just [symlinks](https://en.wikipedia.org/wiki/Symbolic_link) to a
   file named `libffi.so.6.0.4`.

3. Depending on what version you got the `ImportError` regarding, set up a
   symlink of that name to the `libffi.so*` file that you actually have
   installed. For me, I ran `ln -s libffi.so.6.0.4 libffi.so.7` to fix the
   problem.

Once you've set up the symlink, you can restart the Jupyter notebook and you
should be able to log in. (Obvious caveat:
[libffi](https://en.wikipedia.org/wiki/Libffi) is almost as old as I am and I
don't know how it works)

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

### Additional analyses that use the strainFlye pipeline code directly (`sf-analyses`)

**Why does this folder exist?**
The [strainFlye pipeline code](https://github.com/fedarko/strainFlye) includes some
new functionality that we implemented only recently:
for example, strainFlye now supports transversion decoy contexts for FDR
estimation, and the ability to compute p-values for the longest gap between
mutations in a contig.

To avoid mixing the old analysis notebooks (that use _ad hoc_ code to do
various tasks) with newer analyses (that use the actual strainFlye pipeline
code), these "later" analyses are stored in the `sf-analyses` folder.

**What's in this folder?**
The `sf-analyses/chicken/` folder contains analyses of the chicken gut dataset;
the `sf-analyses/sheep/` folder contains analyses of the sheep gut dataset.

**Analysis "starting points."**
Please that the chicken gut analyses begin from the reads and contigs (so, we
first perform alignment using `strainFlye align`).

However, some of the sheep gut analyses here make use of the alignment we
produced in `main-workflow/RUN-ME.sh` above. This 1) avoids rerunning
minimap2 on this massive dataset, and 2) allows us to be consistent with the
rest of these analyses, by using the exact same alignment.

#### Installing strainFlye to run the analyses in `sf-analyses`

Please see the installation instructions in
[strainFlye's README](https://github.com/fedarko/strainFlye).
It's probably easiest to just create a new conda environment and install strainFlye
into that (this is shown in the current strainFlye installation instructions, as of
writing).

If you create a new conda environment to run strainFlye, you will need to
install `jupyter` if you want to run any of the Jupyter notebooks in
`sf-analyses`. As of writing, `jupyter` is included in the main `sheepgut`
environment we set up above in this README -- but it is not included by default
in the strainFlye environment.

(If you got an "internal server error" when
using Jupyter to run the other notebooks, you may get it again in this new
environment -- see the
_Troubleshooting an "internal server error" when running Jupyter notebooks_
section above for some advice.)

### Creating other figures in the paper

The above instructions create *most* of the figures shown in our paper,
although there are a few exceptions.

- strainFlye pipeline figure: created using LibreOffice Draw, GIMP, and Graphviz

- Assembly graph figures: created using
  [MetagenomeScope](https://github.com/marbl/MetagenomeScope) and
  [Bandage](https://github.com/rrwick/Bandage)

- Link graph figure: created using Graphviz (sfdp)
