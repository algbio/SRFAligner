# vcf-to-hapl-to-efg experiment
Pipeline to build an indexable Elastic Founder Graph from a VCF file plus reference. After checking out the 'Prerequisites' and 'Datasets and obtaining the input data' sections of this document, you can build the chromosome 22 iEFG with command
```console
/usr/bin/time ./sample-and-build-efg-heuristic.sh -f chr22_uppercase.fasta -v 1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz -c chr22 -s 2504 -M 8 -t 64
```
This requires ~600 GB of disk space and ~600 GB of RAM. If you want to generate an iEFG from fewer haplotypes, for example 20, run command
```console
/usr/bin/time ./sample-and-build-efg-heuristic.sh -f chr22_uppercase.fasta -v 1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz -c chr22 -s 20 -t 64
```
Finally, the iEFG can be stripped of its paths with command
```console
grep -v "^P" output/efg-unsimplified.gfa > chr22_iEFG.gfa
```

## Prerequisites
The pipeline expects [`bcftools`](https://www.htslib.org/download/) and [`vcf2multialign`](https://github.com/tsnorri/vcf2multialign) to be found in the search path variable `PATH`, and expects `founderblockgraph` to be in folder `tools/founderblockgraph` from the root of this repository. You can get and compile `founderblockgraph` with
```console
git submodule update --init --recursive ../../tools/founderblockgraphs
make -C ../../tools/founderblockgraphs
```

In case you obtain `bcftools` and `vcf2multialign` in a different way, modify lines 11-13 of `sample-and-build-efg-heuristic.sh` accordingly. To manipulate FASTQ files in the next section, we also use [`seqtk`](https://github.com/lh3/seqtk).

## Datasets and obtaining the input data
We use the [phased T2T 1KGP panel](https://zenodo.org/records/7612953) (Version 1.0) by Joseph Lalli, based on [T2T-CHM13v2.0](https://github.com/marbl/CHM13). We can easily obtain chromosome 22 with
```console
wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
seqtk subseq chm13v2.0.fa.gz <(echo "chr22") | seqtk seq -U - > chr22_uppercase.fasta
```
and obtain the chr22 variations using commands
```console
wget "https://zenodo.org/records/7612953/files/phased_T2T_panel.tar"
tar -xvf phased_T2T_panel.tar phased_T2T_panel/1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz --strip-components 1
tar -xvf phased_T2T_panel.tar phased_T2T_panel/1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz.tbi --strip-components 1
```

## iEFG validation
You can validate that the iEFG was correctly built from the MSA rows with `efg-locate`:
```console
../../tools/efg-locate/efg-locate -t 64 --overwrite \
    output/efg-unsimplified.gfa \
    <(awk '{if (substr($0, 0, 1) == ">") {print} else {gsub(/-/, "", $0); print $0}}' \
        output/sampled_haplotypes.a2m) \
    /dev/null
```

## Versions of the software used
| Tool              |          Version |
| ----------------- | ---------------- |
| founderblockgraph | 439ef67 (GitHub) |
| vcf2multialign    |    1.2.2 23f3f42 |
| seqtk             |   1.4-r130-dirty |
| bcftools          |             1.20 |

## todo
- [] gzip some of the intermediate files
