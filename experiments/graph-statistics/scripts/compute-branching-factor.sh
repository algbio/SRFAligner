#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) dag.gfa" ; exit 1
fi

tempdir=$(mktemp -d)
>&2 echo "generated temp dir $tempdir (should be removed automatically on exit)..."
trap '{ rm -rf -- "$tempdir"; }' EXIT

# sort nodes and edges in topological order
grep "^L" $1 | awk '{print $2,$4}' > "$tempdir/edges"
tsort "$tempdir/edges" > "$tempdir/topological_nodes"

awk '(NR == FNR) { inneighbors[$2]=inneighbors[$2] FS $1; next } { \
	n=split(inneighbors[$1],array); \
	for (i=1;i<=n;++i) \
		print array[i],$1}' \
	"$tempdir/edges" "$tempdir/topological_nodes" > "$tempdir/topological_edges"

# find all branching nodes
awk '{outdegree[$1]+=1;} \
	END \
	{ \
		for (key in outdegree) { \
			if (outdegree[key] > 1) {print key} \
		}; \
	};' "$tempdir/edges" > "$tempdir/branching"

# use branching nodes and topological edges to copute the branching factor
gawk --bignum '(NR == FNR) { branching[$1]=1; next } { \
	if (bfactor[$1] + branching[$1] > bfactor[$2]) \
		{bfactor[$2]=bfactor[$1] + branching[$1]} \
	} END { \
		max=0 ; \
		for (key in bfactor) { \
			if (bfactor[key] > max) { \
				max=bfactor[key] \
			} \
		} ;
		print max
	}' "$tempdir/branching" "$tempdir/topological_edges"
