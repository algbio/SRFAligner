# Graph statistics
Analisys of iEFGs or DAGs in GFA format containing only forward `L` links. Command `./run.sh graph.gfa` computes:

- the number of nodes, edges, and base pairs of the graph
- the ⌈N50⌉ metric, that is, the smallest k such that ≥50% of the bases are covered by segments of length ≤k
- the length of the longest node
- the maximum number H of nodes in a block,  if `graph.gfa` is an iEFG
- the width (size of smallest path set covering the nodes, using `GraphChainer`)
- the number of branching nodes (outdegree ≥ 2), choices (sum of outdegrees ≥ 2), the branching factor (maximum number of branching nodes in any path), and the number of maximal paths

## Prerequisites
The scripts used depend on [`octave-cli`](https://octave.org/), `gawk`, `awk`, and `GraphChainer`: the last one is expected to be found in folder `tools/GraphChainer/bin` from the root of this repository, and can be obtained with command
```console
git submodule update --init --recursive ../../tools/GraphChainer
```
and by following the compilation instructions in its README.

## Limitations
The computation of the ⌈N50⌉ is quite slow and inefficient. The computation of the number of maximal paths is also not efficient and keeps in memory many large numbers.
