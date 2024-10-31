#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) efg.gfa" ; exit 1
fi

head -n 3 $1 | grep "^B" | tr "\t" "\n" | tail -n +2 | sort -nr | head -n 1
