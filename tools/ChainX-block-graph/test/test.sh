#!/bin/bash

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
outputfolder=$thisfolder/output-$(date -Iminutes)
logfile=$outputfolder/log
chainxblockgraph=$thisfolder/../chainx-block-graph

mkdir $outputfolder
echo -n > $logfile

for testfile in $thisfolder/input/anchors-1.gaf
do
	echo "$testfile : " >> $logfile
	basename=$(basename $testfile)
	outfileg=$outputfolder/${basename%.*}-global.gaf
	correctg=$thisfolder/correctoutput/${basename%.*}-global.gaf
	$chainxblockgraph --unsorted-input --global $thisfolder/input/graph.gfa $testfile $outfileg \
		>> $logfile 2>> $logfile
	diff <(sort $outfileg) <(sort $correctg) > /dev/null 2>/dev/null
	exitcode=$? ; if [ $exitcode -ne 0 ] ; then
		echo "Test failed on file $testfile!" | tee -a $logfile
		exit 1
	fi

	outfilesg=$outputfolder/${basename%.*}-semi-global.gaf
	correctsg=$thisfolder/correctoutput/${basename%.*}-semi-global.gaf

	$chainxblockgraph --unsorted-input --semi-global $thisfolder/input/graph.gfa $testfile $outfilesg \
		>> $logfile 2>> $logfile

	diff <(sort $outfilesg) <(sort $correctsg) > /dev/null 2>/dev/null
	exitcode=$? ; if [ $exitcode -ne 0 ] ; then
		echo "Test failed on file $testfile!" | tee -a $logfile
		exit 1
	fi
done
