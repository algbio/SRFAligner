#!/bin/bash
set -eo pipefail
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) graph.gfa" ; exit 1
fi

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
graphchainer=$thisfolder/../../../tools/GraphChainer/bin/GraphChainer

$graphchainer --graph-statistics -g $1 -f fakereads.fastq -a fakealns.gaf 2>&1 | tail -n 1 | cut -d' ' -f8
