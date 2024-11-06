# vcf-to-hapl-to-efg experiment
Pipeline to build an indexable Elastic Founder Graph from a VCF file plus reference. After checking out the 'Prerequisites' and 'Datasets and obtaining the input data' sections, you can build the chromosome 22 iEFG with command
```
/usr/bin/time ./sample-and-build-efg-heuristic.sh -f chr22_uppercase.fasta -v 1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz -c chr22 -s 2504 -M 250 -t 32
```
This requires ~X of disk space and ~Y GB of RAM. If you want to generate an iEFG from fewer haplotypes, for example 20, run
```
/usr/bin/time ./sample-and-build-efg-heuristic.sh -f chr22_uppercase.fasta -v 1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz -c chr22 -s 20 -M 41 -t 32
```

## Prerequisites
The pipeline expects [`bcftools`](https://www.htslib.org/download/) and [`vcf2multialign`](https://github.com/tsnorri/vcf2multialign) to be found in the search path variable `PATH`, and expects `founderblockgraph` to be in folder `tools/founderblockgraph` from the root of this repository. In case you obtain the tools in a different way, modify lines 11-33 of `sample-and-build-efg-heuristic.sh` accordingly. To manipulate FASTQ files in the next section, we also use [`seqtk`](https://github.com/lh3/seqtk).

## Datasets and obtaining the input data
We use the [phased T2T 1KGP panel](https://zenodo.org/records/7612953) (Version 1.0) by Joseph Lalli, based on [T2T-CHM13v2.0](https://github.com/marbl/CHM13). Assuming `seqtk` is installed, we can easily obtain chromosome 22 as follows:
```
wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
seqtk subseq chm13v2.0.fa.gz <(echo "chr22") | seqtk seq -U - > chr22_uppercase.fasta
```
and obtain the chr22 variations as follow:
```
wget "https://zenodo.org/records/7612953/files/phased_T2T_panel.tar"
tar -xvf phased_T2T_panel.tar phased_T2T_panel/1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz --strip-components 1
tar -xvf phased_T2T_panel.tar phased_T2T_panel/1KGP.CHM13v2.0.chr22.recalibrated.snp_indel.pass.phased.vcf.gz.tbi --strip-components 1
```

## Versions of the software used
| Tool              |          Version |
| ----------------- | ---------------- |
| founderblockgraph | 3f4133e (GitHub) |
| vcf2multialign    |    1.2.2 23f3f42 |
| seqtk             |   1.4-r130-dirty |
| bcftools          |             1.20 |

## todo
- [] gzip some of the intermediate files
