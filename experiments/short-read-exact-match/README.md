# Exact match of short reads on the chr22 iEFG
We compare the short-read exact matching solution of `efg-locate` on the chromosome 22 iEFG built with the pipeline at `experiments/vcf-to-hapl-to-efg` to that of: `bwa`, on the T2T-CHM13 linear reference for chromosome 22; and `vg map`, on the (pruned) chromosome 22 graph built from the same VCF as the iEFG. After checking out the *Prerequisites* and *Datasets* sections, run the script `runexp.sh` (requires ~150G of disk space for the results) and check `output/runexp_log.txt` for the results.

## Prerequisites
Script `runexp.sh` expects `efg-locate`, `bwa`, `seqtk`, and `vg` to be located in folders `tools/efg-locate`, `tools/bwa`, `tools/seqtk`, and `tools/vg/bin` from the root of this repository. You can download and compile them with the following commands (executed from this folder):
```console
make -C ../../tools/efg-locate
git submodule update --init ../../tools/{bwa,seqtk}
mkdir --parents ../../tools/vg/bin && wget https://github.com/vgteam/vg/releases/download/v1.67.0/vg --output-document=../../tools/vg/bin
make -C ../../tools/bwa
make -C ../../tools/seqtk
```

## Datasets
Obtain the graphs (544M) with command
```console
wget https://zenodo.org/records/15189487/files/chr22_iEFG.gfa.gz?download=1 --output-document=input/chr22_iEFG.gfa.gz && gunzip input/chr22_iEFG.gfa.gz
wget https://zenodo.org/records/15189487/files/chr22_vg.gfa.gz?download=1 --output-document=input/chr22_vg.gfa.gz && gunzip input/chr22_vg.gfa.gz
```
the reads (24GB) with
```console
wget 'https://cs.helsinki.fi/group/gsa/panvc-founders/scalability-experiment/reads/ERR1025645_sample05_1.fq.gz' --output-document=input/ERR1025645_sample05_1.fq.gz
```
and the chromosome 22 reference (1GB) with
```console
wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" --output-document=input/chm13v2.0.fa.gz
seqtk subseq input/chm13v2.0.fa.gz <(echo "chr22") | seqtk seq -U - > input/chr22_uppercase.fasta
```
