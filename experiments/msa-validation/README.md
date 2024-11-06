# Evaluation of sequence-to-graph aligners on simulated reads
Script to compare the output of `vcf2multialign -H` on a chromosome of the [T2T-CHM13v2.0](https://github.com/marbl/CHM13) + [Phased T2T 1KGP panel](https://zenodo.org/records/7612953#.Y-8VD3bMJPY) dataset with that of [`bcftools consensus`](https://samtools.github.io/bcftools/howtos/consensus-sequence.html). Specifically, each sequence of the MSA generated with `vcf2multialign` is stripped of gaps and compared to the corresponding output of `bcftools consensus` using [`edlib-aligner`](https://github.com/Martinsos/edlib). Each Levenshtein distance (NW) is collected and the min, max, and average distance are computed.

## Prerequisites
The script depends on `awk` and `bcftools`. It expects to find executable `edlib-aligner` and `seqtk` in folders `tools/edlib/meson-build` and `tools/seqtk/seqtk` from the root of the repository, which can be compiled as follows:
```console
git submodule update --init ../../tools/edlib
make -C ../../tools/edlib
git submodule update --init ../../tools/seqtk
make -C ../../tools/seqtk
```

## Usage
Usage is as follows:
```console
./validate.sh MSA.fasta reference.fasta variations.vcf.gz threads output_stats.txt
```
where `MSA.fasta` is the MSA computed by `vcf2multialign -H`, `reference.fasta` contains the reference chromosome, `variations.vcf.gz` contains the variations to the chromosome, threads is a positive number of threads, and `output_stats.txt` is the desired output file. For the chromosome 22 built in experiment `vcf-to-hapl-to-efg`, the command to run is
```console
./validate.sh ../vcf-to-hapl-to-efg/output/sampled_haplotypes.a2m ../vcf-to-hapl-to-efg/chr22_uppercase.fasta ../vcf-to-hapl-to-efg/1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz 64 output_stats.txt
```
