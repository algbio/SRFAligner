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
# following GraphChainer compile instructions
cd ../../tools/GraphChainer
conda env create -f CondaEnvironment.yml
conda activate GraphChainer
make bin/GraphChainer

cd ../../experiments/aligner-evaluation
conda activate aligner-evaluation
```
The scripts `runexp.sh` expect `/usr/bin/time`, `openssl`, and `awk` to be installed (and the `aligner-evaluation` conda environment to be active).

## Datasets
| Graph                       | Construction                           | Input |
|-----------------------------|----------------------------------------|-------|
| chr22\_iEFG                 | `vcf2multialign` + `founderblockgraph` | [T2T-CHM13v2.0](https://github.com/marbl/CHM13) + [Phased T2T 1KGP panel](https://zenodo.org/records/7612953#.Y-8VD3bMJPY)
| chr22\_vg                   | `vg`                                   | same as above |
| chr22_vg_msa                | `vcf2multialign` + `vg -M`             | same as above |
| covid19_100_iEFG_simplified | `founderblockgraph` + `efg-simplify`   | [MSA from efg-mems experiments](https://github.com/algbio/efg-mems)

All experiments except for `vg-comparison` and `mems` use the chromosome 22 iEFG built with the `vcf-to-hapl-to-efg` pipeline as file `input/chr22_iEFG.gfa`. You can get the graph from [zenodo](https://doi.org/10.5281/zenodo.14012882) as follows:
```
wget https://zenodo.org/records/14012882/files/chr22_iEFG.gfa.gz?download=1 --output-document=input/chr22_iEFG.gfa.gz
gunzip input/chr22_iEFG.gfa
```
The `mems` experiment uses a SaRS-CoV-2 MSA of 100 strains (NCBI accession numbers in `input/covid19_100_acc.txt`) aligned with [ViralMSA](https://github.com/niemasd/ViralMSA):
```
wget https://zenodo.org/records/14012882/files/covid19_100_iEFG_simplified.gfa.gz?download=1 --output-document=input/covid19_100_iEFG_simplified.gfa.gz
gunzip input/covid19_100_iEFG_simplified.gfa.gz
```
Finally, the `vg-comparison` experiment additionally uses the `vg` graphs built from the 1KGP phased VCF file or the MSA obtained with `vcf2multialign` from the same data:
```
wget https://zenodo.org/records/14012882/files/chr22_vg.gfa.gz?download=1 --output-document=input/chr22_vg.gfa.gz
gunzip input/chr22_vg.gfa.gz
wget https://zenodo.org/records/14012882/files/chr22_vg_msa.gfa.gz?download=1 --output-document=input/chr22_vg_msa.gfa.gz
gunzip input/chr22_vg_msa.gfa.gz
```
