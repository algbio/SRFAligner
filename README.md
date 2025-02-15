# seed-chain-extend on indexable Elastic Founder Graphs
`SRFAligner` and `SRFChainer` are long-read alignment tools based on indexable Elastic Founder Graphs (iEFGs). iEFGs can be obtained from FASTA multiple sequence alignments using [`founderblockgraph`](https://github.com/algbio/founderblockgraphs), or from a VCF file with the pipeline implemented in `experiments/vcf-to-hapl-to-efg` using `founderblockgraph` and [`vcf2multialign`](https://github.com/tsnorri/vcf2multialign/). The graphs used in the experiments can be found at [doi.org/10.5281/zenodo.14012881](https://doi.org/10.5281/zenodo.14012881).

![Workflow to build iEFGs from a VCF file and to perform seed-chain-extend alignment](docs/workflow.png)

## getting started
`SRFAligner` and `SRFChainer` depend on `efg-locate`, `chainx-block-graph`, and `GraphAligner`. To clone this repository and compile the first two:
```
git clone https://github.com/algbio/SRFAligner && cd SRFAligner
git submodule update --init tools/{sdsl-lite-v3,concurrentqueue}
make
```
and `GraphAligner`'s executable is expected to be found in folder `tools/GraphAligner/bin`, so you can run command `git submodule update --init --recursive tools/GraphAligner` and follow its [compilation instructions](https://github.com/maickrau/GraphAligner?tab=readme-ov-file#compilation). If `GraphAligner` is already installed in your system, you can just modify the relative line in `SRFAligner` and `SRFChainer`:
```
sed --in-place '7s/.*/graphaligner=GraphAligner/' SRFAligner
sed --in-place '7s/.*/graphaligner=GraphAligner/' SRFChainer
sed --in-place '7s/.*/graphaligner=GraphAligner/' efg-memsAligner
```
Finally, to use MEM seeds computed by `efg-mems`, `efg-memsAligner` expects `efg-mems`'s executable to be in `tools/efg-mems/efg-mems`:
```
git submodule update --init --recursive tools/efg-mems
cd tools/efg-mems/sdsl-lite
./install.sh .
cd ..
cmake .
```
