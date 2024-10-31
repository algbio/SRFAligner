#!/bin/bash
set -eo pipefail
if [[ $# -lt 2 ]] || [[ $# -gt 2 ]]
then
	>&2 echo "Usage: $(basename $0) graph.gfa [0-100]" ; exit 1
fi

# Configuration
quantiles="$(echo $2 / 100 | bc -l)"

lengths=$(mktemp)
>&2 echo "generated tmp file $lengths (should be removed automatically on exit)..."
trap '{ rm -f -- "$lengths"; }' EXIT

grep "^S" $1 | awk '{for (i = 1 ; i <= length($3) ; i++) {print length($3)}}' | sort -n > $lengths

octaveout=$(octave-cli --eval "format long; x = dlmread('$lengths'); q = quantile(x, [0.00 $quantiles]); N = ceil(q([2:length(q)])); disp(N)" | tr -s " " "\t" | cut -f2-)

echo -e "$octaveout"
