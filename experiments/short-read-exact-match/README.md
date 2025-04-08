# Exact match of short reads on the chr22 iEFG
We compare the short-read exact matching solution of `efg-locate` on the chromosome 22 iEFG built with the pipeline at `experiments/vcf-to-hapl-to-efg` to that of `bwa` on the T2T-CHM13 linear reference for chromosome 22. After checking out the *Prerequisites* and *Datasets* sections, run the script `runexp.sh` (requires 42G of disk space) and check `output/runexp_log.txt` for the results.

## Prerequisites
Script `runexp.sh` expects `efg-locate`, `bwa`, and `seqtk` to be located in folder `tools/efg-locate`, `tools/bwa`, and `tools/seqtk` from the root of this repository. You can download and compile them with the following commands (executed from this folder):
```console
make -C ../../tools/efg-locate
git submodule update --init ../../tools/{bwa,seqtk}
make -C ../../tools/bwa
make -C ../../tools/seqtk
```

## Datasets
Obtain the graph (262M) with command
```console
wget https://zenodo.org/records/15112649/files/chr22_iEFG.gfa.gz?download=1 --output-document=input/chr22_iEFG.gfa.gz && gunzip input/chr22_iEFG.gfa.gz
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
