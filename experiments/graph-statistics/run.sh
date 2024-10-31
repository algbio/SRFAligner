#!/bin/bash
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]
then
	>&2 echo "Usage: $(basename $0) graph.gfa" ; exit 1
fi

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script

echo -e "nodes:\t$($thisfolder/scripts/compute-nodes.sh $1)"
echo -e "edges:\t$($thisfolder/scripts/compute-edges.sh $1)"
echo -e "bases:\t$($thisfolder/scripts/compute-bps.sh $1)"
echo -e "N50:\t$($thisfolder/scripts/compute-N.sh $1 50 2>> /dev/null)" 
echo -e "longest node:\t$($thisfolder/scripts/compute-longest-node.sh $1)"
echo -e "H:\t$($thisfolder/scripts/compute-efg-H.sh $1)"
echo -e "width:\t$($thisfolder/scripts/compute-width.sh $1)"
echo -e "branching nodes:\t$($thisfolder/scripts/compute-branching-nodes.sh $1)"
echo -e "choices:\t$($thisfolder/scripts/compute-choices.sh $1)"
echo -e "branching factor:\t$($thisfolder/scripts/compute-branching-factor.sh $1 2>> /dev/null)"
echo -e "paths:\t$($thisfolder/scripts/compute-paths.sh $1 2>> /dev/null)"
