#!/bin/bash
# remember to run `conda activate aligner-evaluation`
set -e
set -o pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

# executable's absolute paths/commands (make sure they work!)
BADREAD=$thisfolder/../../../tools/Badread/badread-runner.py
GRAPHCHAINER=$thisfolder/../../../tools/GraphChainer/bin/GraphChainer
graphaligner=$thisfolder/../../../tools/GraphAligner/bin/GraphAligner
usrbintime=/usr/bin/time

# params
inputgraphs=($thisfolder/../input/t2t-1KGP-phased-hapl-complete-250-heur-unsimplified.gfa $thisfolder/../input/vg-t2t-1KGP-phased-hapl-complete.gfa $thisfolder/../input/vgmsa-t2t-1KGP-phased-hapl-complete.gfa)
graphnames=(iefg vg vgmsa)
coverage=30
threads=64

# 0. setup
mkdir output
echo -n > output/runexp_log.txt

# 1. simulate path and reads
for ((g = 0 ; g < ${#inputgraphs[@]} ; g++))
do
	# TODO maybe parallelize this?
	inputgraph="${inputgraphs[$g]}"
	graphname="${graphnames[$g]}"
	python3 ../scripts/generate_sim_reads.py \
		--graph $inputgraph \
		--fastq output/sim_reads_$graphname.fastq \
		--seed $(openssl enc -aes-256-ctr -pass pass:"semirepeatfree" -nosalt </dev/zero 2>/dev/null | shuf -i0-4294967295 -n 1 --random-source=/dev/stdin) \
		--path output/sim_reads_path_$graphname.nodes \
		--fasta output/sim_reads_path_$graphname.fasta \
		--coverage $coverage \
		2>> output/runexp_log.txt >> output/runexp_log.txt
done

# 3. run the aligners on each dataset
for ((g1 = 0 ; g1 < ${#inputgraphs[@]} ; g1++))
do
	for ((g2 = 0 ; g2 < ${#inputgraphs[@]} ; g2++))
	do
		inputgraph="${inputgraphs[$g1]}"
		graphname="${graphnames[$g1]}"
		datasetname="${graphnames[$g2]}"
		reads=output/sim_reads_$datasetname.fastq
		$usrbintime $graphaligner \
			-t $threads \
			-x vg \
			-g $inputgraph \
			-f $reads \
			-a output/${graphname}_graph_${datasetname}_reads_alignments.gaf \
			2>> output/runexp_log.txt >> output/runexp_log.txt

		# pick first alignment
		awk '{if (found[$1] == "1") \
		          {} \
		      else
		          {found[$1]="1"; print}}' \
		output/${graphname}_graph_${datasetname}_reads_alignments.gaf \
		> output/best_${graphname}_graph_${datasetname}_reads_alignments.gaf
	done
done

# 4. validate and plot results
for ((g1 = 0 ; g1 < ${#inputgraphs[@]} ; g1++))
do
	for ((g2 = 0 ; g2 < ${#inputgraphs[@]} ; g2++))
	do
		inputgraph="${inputgraphs[$g1]}"
		graphname="${graphnames[$g1]}"
		if [ $g1 -eq $g2 ]
		then
			# we have ground truth
			reads=output/sim_reads_$graphname.fastq
			path=output/sim_reads_path_$graphname.nodes
			fasta=output/sim_reads_path_$graphname.fasta
			alignments=output/best_${graphname}_graph_${graphname}_reads_alignments.gaf
			python3 ../scripts/compute_summary.py \
				-t 3 \
				--graph $inputgraph \
				--fastq $reads \
				--path $path \
				--fasta $fasta \
				--alignments $alignments \
				--metrics output/metrics_${graphname}_graph_${graphname}_reads.mts &
		else
			# we do not have ground truth
			datasetname="${graphnames[$g2]}"
			reads=output/sim_reads_$datasetname.fastq
			alignments=output/best_${graphname}_graph_${datasetname}_reads_alignments.gaf
			python3 ../scripts/compute_summary.py \
				-t 3 \
				--graph $inputgraph \
				--fastq $reads \
				--alignments $alignments \
				--metrics output/metrics_${graphname}_graph_${datasetname}_reads.mts &
		fi
	done
done
wait $(jobs -p)

for ((g = 0 ; g < ${#inputgraphs[@]} ; g++))
do
	graphname="${graphnames[$g]}"
	python3 ../scripts/compute_metrics.py \
		--output-name output/results_${graphname}_graph \
		--summaries output/metrics_${graphname}_graph_*.mts \
		--summaries-names $(basename -s ".mts" -a output/metrics_${graphname}_graph_*.mts | sed 's/metrics_//')
done
