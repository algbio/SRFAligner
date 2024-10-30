# seed-chain-extend on indexable Elastic Founder Graphs
`SRFAligner` and `SRFChainer` are long-read alignment tools based on indexable Elastic Founder Graphs (iEFGs). iEFGs can be obtained from FASTA multiple sequence alignments using [`founderblockgraph`](https://github.com/algbio/founderblockgraphs), or from a VCF file with the pipeline implemented in `experiments/vcf-to-hapl-to-efg` using `founderblockgraph` and [`vcf2multialign`](https://github.com/tsnorri/vcf2multialign/).

![Workflow to build iEFGs from a VCF file and to perform seed-chain-extend alignment](docs/workflow.png)

## getting started
`SRFAligner` and `SRFChainer` depend on `efg-locate`, `chainx-block-graph`, and `GraphAligner`. To clone this repository and compile the first two:
```
git clone https://github.com/algbio/SRFAligner && cd SRFAligner
git submodule update --init tools/sdsl-lite-v3
make
```
and `GraphAligner`'s executable is expected to be found in folder `tools/GraphAligner/bin`, so you can execute `git submodule update --init --recursive tools/GraphAligner` and follow its compilation instructions. If `GraphAligner` is already installed in your system, you can just modify the relative line in `SRFAligner` and `SRFChainer`:
```
sed --in-place '7s/.*/graphaligner=GraphAligner/' SRFAligner
sed --in-place '7s/.*/graphaligner=GraphAligner/' SRFChainer
```
