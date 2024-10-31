#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) graph.gfa" ; exit 1
fi

grep "^S" $1 | awk '{if (length($3) > max) {max = length($3)}} END {print max}'
