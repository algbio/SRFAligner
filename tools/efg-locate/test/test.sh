#!/bin/bash
# barebones testing pipeline that matches the test files with the correct
# outputs as specified by the following arrays
locate=("tcs_fig_5.gfa tcs_fig_5_edge.fasta        tcs_fig_5_edge.gfa"
	"tcs_fig_5.gfa tcs_fig_5_three_nodes.fasta tcs_fig_5_three_nodes.gfa"
	"tcs_fig_5.gfa tcs_fig_5_four_nodes.fasta  tcs_fig_5_four_nodes.gfa"
	"indels.gfa    indels_five_nodes.fasta     indels_five_nodes.gfa")

approximate=("tcs_fig_5.gfa tcs_fig_5_approximate.fasta tcs_fig_5_approximate.fasta")

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
outputfolder=$thisfolder/output-$(date -Iminutes)
logfile=$outputfolder/log
efglocate=$thisfolder/../efg-locate

mkdir $outputfolder
echo -n > $logfile

for testfile in "${locate[@]}"
do
	graph=$thisfolder/inputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f1)
	patterns=$thisfolder/inputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f2)
	correct=$thisfolder/outputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f3)

	patternsbasename=$(basename $patterns)
	output=$outputfolder/${patternsbasename%.*}.gfa

	echo "$efglocate $graph $patterns $output" >> $logfile
	$efglocate $graph $patterns $output >> $logfile 2>> $logfile
	diff $output $correct > /dev/null 2>/dev/null

	exitcode=$? ; if [ $exitcode -ne 0 ] ; then
		echo "Test failed for files $graph $patterns $correct!" | tee -a $logfile
		exit 1
	fi
done

for testfile in "${approximate[@]}"
do
	graph=$thisfolder/inputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f1)
	patterns=$thisfolder/inputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f2)
	correct=$thisfolder/outputs/$(echo "$testfile" | tr -s " " | cut -d' ' -f3)

	patternsbasename=$(basename $patterns)
	output=$outputfolder/${patternsbasename%.*}.gaf

	echo "$efglocate $graph $patterns $output" >> $logfile
	$efglocate --approximate $graph $patterns $output >> $logfile 2>> $logfile
	diff $output $correct > /dev/null 2>/dev/null

	exitcode=$? ; if [ $exitcode -ne 0 ] ; then
		echo "Test failed for files $graph $patterns $correct!" | tee -a $logfile
		exit 1
	fi
done
