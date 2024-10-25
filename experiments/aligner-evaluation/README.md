# Evaluation of sequence-to-graph aligners on simulated reads
Based on the [GraphChainer evaluation pipeline](https://github.com/algbio/GraphChainer-scripts).

## Prerequisites
Install [miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Then, set up the environment `aligner-evaluation` defined in `environment.yml`:
```
conda env create -f environment.yml
```
The pipeline expects to find `GraphChainer` and `badread_runner.py` in folders `tools/GraphChainer/bin` and `tools/BadRead` from the root of the repository, so download them and compile `GraphChainer` if you have not.
```
git submodule update --init --recursive ../../tools/Badread
git submodule update --init --recursive ../../tools/GraphChainer
# follow GraphChainer building instructions
```
The scripts `runexp.sh` expect `/usr/bin/time`, `openssl`, and `awk` to be installed.

## Datasets
todo

## todo
- [] docs
- [] test GraphAligner on vg graphs
- [] test other aligners
