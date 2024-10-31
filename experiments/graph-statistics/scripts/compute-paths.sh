#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) connected-dag.gfa" ; exit 1
fi

tempdir=$(mktemp -d)
>&2 echo "generated temp dir $tempdir (should be removed automatically on exit)..."
trap '{ rm -rf -- "$tempdir"; }' EXIT

# compute topological nodes and edges
grep "^L" $1 | awk '{print $2,$4}' > "$tempdir/edges"
tsort "$tempdir/edges" > "$tempdir/topological_nodes"

awk '(NR == FNR) { inneighbors[$2]=inneighbors[$2] FS $1; next } { \
	n=split(inneighbors[$1],array); \
	for (i=1;i<=n;++i) \
		print array[i],$1}' \
	"$tempdir/edges" "$tempdir/topological_nodes" > "$tempdir/topological_edges"

# find sinks in the graph
awk '{nodes[$2]=1; nonsinks[$1]=1} \
	END \
	{ \
		for (key in nodes) { \
			if (nonsinks[key] != 1) {print key} \
		}; \
	};' "$tempdir/edges" > "$tempdir/sinks"

# use edges and sinks to compute paths to each node
gawk --bignum '(NR == FNR) { sinks[$1]=1; next } { \
	if (paths[$1] > 0) \
		{paths[$2]+=paths[$1]} \
	else \
		{paths[$2]+=1}}
	END { \
		result=0 ; \
		for (key in sinks) { \
			result+=paths[key]} ; \
		printf "%e", result ; \
		print ""\
	}' "$tempdir/sinks" "$tempdir/topological_edges"
