#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) graph.gfa" ; exit 1
fi

grep "^L" $1 | \
	awk '{outdegree[$2]+=1;} \
	END \
	{ \
		for (key in outdegree) { \
			if (outdegree[key] > 1) {result+=1} \
		}; \
		print result \
	};'
