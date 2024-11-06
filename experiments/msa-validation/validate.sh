#!/bin/bash
set -eo pipefail
if [[ $# -lt 5 ]] || [[ $# -gt 5 ]]
then
	>&2 echo "Usage: $(basename $0) MSA.fasta reference.fasta variations.vcf.gz threads output_stats.txt" ; exit 1
fi

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
export bcftools=bcftools
export edlib=../../tools/edlib/meson-build/edlib-aligner
export seqtk=../../tools/seqtk/seqtk

msa=$1
export ref=$2
export vcf=$3
threads=$4
outputfile=$5

if [ -f "$outputfile" ]
then
	>&2 echo "$outputfile already exists!"; exit 1
fi

function align_destroy () { # $1 is temp file containing gapless fasta sequence, to be destroyed
	header=$(head -n 1 $1)
	sample=$(echo ${header:1} | cut -d'-' -f1)
	haplot=$(echo ${header:1} | cut -d'-' -f2)
	echo -en "$sample-$haplot\t"
	$edlib \
		$1 \
		<($bcftools consensus -f $ref -s $sample -H $haplot $vcf 2>> /dev/null | $seqtk seq -U /dev/stdin) \
		| grep "^#0" | cut -d' ' -f2
	rm $1
}
# make align_destroy visible to GNU parallel
export SHELL=$(type -p bash)
export -f align_destroy

while read -u 3 header
do
	while [[ $(find . -type f -name "tmp*" | wc -l) -gt $threads ]]
	do
		sleep 1s
	done
	read -u 3 entry
	sample=$(echo ${header:1} | cut -d'-' -f1)
	haplot=$(echo ${header:1} | cut -d'-' -f2)
	if [ "$sample" == "REF" ] ; then continue ; fi
	gaplessentry=$(echo "$entry" | tr -d "-")

	tempfile=$(mktemp -p $thisfolder)
	echo "$header" > $tempfile
	echo "$gaplessentry" >> $tempfile
	echo "$tempfile"
done 3<$msa | parallel --keep-order --jobs $threads align_destroy > $outputfile

awk 'BEGIN { min = 2^1024 ; max = 0 } \
	{ sum += $2 ; count++ ; \
	  if ($2 < min) min = $2 ; \
	  if ($2 > max) max = $2 } \
	END { print "average:\t"sum / count ; print "    min:\t"min ; print "    max:\t"max }' \
	$outputfile >> $outputfile

tail -n 3 $outputfile

