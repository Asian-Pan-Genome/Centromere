# HORmining

## Overview

HORmining identifies higher-order repeat structures in decomposed alpha-satellite monomers. It provides two complementary methods:

- `HORmon`: graph-based HOR identification.
- `HiCAT`: hierarchical tandem-repeat mining.

This portable revision preserves the original scientific algorithms and focuses on reproducible installation, clearer command-line validation, headless server execution, and a minimal smoke-test dataset.

## Changes in this revision

- Replaced the machine-specific Conda export with a portable, directly pinned `environment.yml`.
- Standardized the environment name as `HORminingEnv`.
- Removed the original absolute Conda `prefix` path.
- Added Graphviz and PyGraphviz, which are required by the HORmon graph-output code.
- Restored the Levenshtein and Matplotlib imports used by `HiCAT_HOR.py` and selected the non-interactive `Agg` backend for servers.
- Added required-argument validation, method choices, an integer distance parameter, and `--version`.
- Added deterministic minimal test data and a smoke-test script for both methods.

## Supported platform

The supplied environment is intended for Linux servers and HPC systems. It does not require a graphical desktop. Conda, Mamba, or Micromamba is required to create the environment.

## Installation

From the pipeline directory:

```bash
conda env create -f environment.yml
conda activate HORminingEnv
```

If Mamba is available, it is generally faster:

```bash
mamba env create -f environment.yml
conda activate HORminingEnv
```

Verify the installation:

```bash
python HORmining.py --version
bedtools --version
dot -V
python -c "import pybedtools, pygraphviz, networkx, numpy, pandas; print('Python dependencies OK')"
```

Expected version output:

```text
HORmining 1.0.0
```

## Minimal server smoke test

The `tests/data` directory contains a small artificial monomer array consisting of repeated `1-2-3` units (36 monomers; 12 repeats). Cluster identifiers must be integers because the summary code normalizes HORs numerically. This dataset is intended only to confirm installation and execution; it is not a biological benchmark.

Run both methods:

```bash
bash tests/run_smoke_test.sh
```

The script creates a timestamped directory under `tests/smoke_outputs`, runs HiCAT and HORmon, verifies their principal output files, and prints the output location. A successful run ends with:

```text
HORmining smoke tests passed
```

For a permanent record of the server test:

```bash
bash tests/run_smoke_test.sh 2>&1 | tee server_smoke_test.log
```

Please retain `server_smoke_test.log` and the generated output directory when reporting test results.

## Input files

### Decomposed alpha-satellite BED

Four tab-separated columns are required:

```text
chromosome  start  end  monomer_identifier
```

The monomer identifier must encode coordinates and strand in one of these forms:

```text
chr1:100-271(+)
sample::chr1:100-271(-)
```

### Monomer cluster table

Two tab-separated columns are required:

```text
monomer_identifier  cluster_identifier
```

Every monomer identifier in the BED file should have a corresponding integer cluster assignment. Use labels such as `1`, `2`, and `3`, not letters such as `A`, `B`, and `C`.

## Usage

### HiCAT method

```bash
python HORmining.py \
  -bed input.asat.bed \
  -c input.clusters.tsv \
  -O output_hicat \
  -m HiCAT
```

Principal summary:

```text
output_hicat/HiCAT/HiCAT.hor.summary.xls
```

### HORmon method

```bash
python HORmining.py \
  -bed input.asat.bed \
  -c input.clusters.tsv \
  -O output_hormon \
  -m HORmon
```

Principal decomposition summary:

```text
output_hormon/graph/graph.HORdecomposition.stat.xls
```

### Optional distance parameter

`-d` specifies the maximum distance in base pairs allowed between adjacent monomers assigned to the same continuous block. The default is `10000`.

```bash
python HORmining.py -bed input.asat.bed -c input.clusters.tsv -O output -m HORmon -d 10000
```

## Recommended server test checklist

Run the following commands and retain their output:

```bash
uname -a
conda --version
conda env create -f environment.yml
conda activate HORminingEnv
python HORmining.py --version
bedtools --version
dot -V
bash tests/run_smoke_test.sh 2>&1 | tee server_smoke_test.log
```

After the smoke test passes, run one representative real dataset with each method if computationally feasible. Compare the revised outputs with outputs from the original pipeline for the same inputs.

## Recording the fully resolved Linux environment

The directly used dependencies are pinned in `environment.yml`. After a successful server installation, record the exact solved package builds for long-term reproducibility:

```bash
conda list --explicit > conda-linux-64.lock.txt
```

Commit `conda-linux-64.lock.txt` with the tested release. It can later recreate the exact Linux environment with:

```bash
conda create -n HORmining-lock --file conda-linux-64.lock.txt
```

## Troubleshooting

### Python 3.9 and PyGraphviz version conflict

PyGraphviz 1.14 cannot be installed with Python 3.9.18. The revised environment pins `pygraphviz=1.11` and `graphviz=9.0.0` together: the conda-forge Linux CPython 3.9 build of PyGraphviz 1.11 (`py39h7ef11f9_2`) requires Graphviz >=9.0.0,<10.0a0. Changing only PyGraphviz while keeping Graphviz 12.2.1 would leave a conflict.

These pins were checked against package metadata; a complete Linux environment solve and execution test remain to be performed on the server. After uploading the corrected files, first check resolution without installing:

```bash
mamba env create -f environment.yml --dry-run
```

If successful, run the installation and smoke test above. Cache modification warnings are separate from the Python/PyGraphviz version conflict shown by the solver.

### `dot` was not found

Confirm that the `HORminingEnv` environment is active and run `dot -V`. Recreate the environment from `environment.yml` if necessary.

### PyGraphviz import or linking error

Install Graphviz and PyGraphviz from the same Conda channel. The supplied environment uses conda-forge for both packages.

### Bedtools or pybedtools error

Run `bedtools --version` and confirm that the BED file is tab-delimited, coordinate-sorted, and contains four columns.

### Server has no display

No display is required. Matplotlib is configured to use the `Agg` backend.

## Release notes

Before public release, add the project license and archive a tested, versioned release. The minimal smoke test establishes software execution only and does not validate biological accuracy.
