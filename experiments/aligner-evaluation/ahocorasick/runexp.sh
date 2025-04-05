#!/bin/bash
# comparison of SRFAligner and efg-ahocorasickAligner
# remember to run `conda activate aligner-evaluation` before execution

set -e
set -o pipefail

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

# executable's absolute paths/commands (make sure they work!)
graphaligner=$thisfolder/../../../tools/GraphAligner/bin/GraphAligner
srfaligner=$thisfolder/../../../SRFAligner
ahocorasickaligner=$thisfolder/../../../efg-ahocorasickAligner
usrbintime=/usr/bin/time

# params
inputgraph=$thisfolder/../input/chr22_iEFG.gfa
coverage=30
threads=1

# 0. setup
mkdir output
echo -n > output/runexp_log.txt
set -a
set +a

# 1. simulate path and reads
# uncomment the following 3 lines and comment the following ones to use the reads in semi-repeat-free
#ln -s $thisfolder/../semi-repeat-free/output/sim_reads.fastq output
#ln -s $thisfolder/../semi-repeat-free/output/sim_reads_path.fasta output
#ln -s $thisfolder/../semi-repeat-free/output/sim_reads_path.nodes output
python3 ../scripts/generate_sim_reads.py \
	--graph $inputgraph \
	--fastq output/sim_reads.fastq \
	--seed $(openssl enc -aes-256-ctr -pass pass:"semirepeatfree" -nosalt </dev/zero 2>/dev/null | shuf -i0-4294967295 -n 1 --random-source=/dev/stdin) \
	--path output/sim_reads_path.nodes \
	--fasta output/sim_reads_path.fasta \
	--coverage $coverage \
	2>> output/runexp_log.txt >> output/runexp_log.txt

# 2. run the aligners
$usrbintime $srfaligner \
	-t $threads \
	-g $inputgraph \
	-m 0 \
	-p \
	-f output/sim_reads.fastq \
	-a output/semi_repeat_free_alignments.gaf \
	-w output \
	2>> output/runexp_log.txt >> output/runexp_log.txt

$usrbintime $ahocorasickaligner \
	-g $inputgraph \
	-f output/sim_reads.fastq \
	-a output/ahocorasick_alignments.gaf \
	-w output \
	2>> output/runexp_log.txt >> output/runexp_log.txt

# 3. pick first alignment
for alignment in "ahocorasick_alignments.gaf" "semi_repeat_free_alignments.gaf"
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/$alignment > output/best_$alignment
done

# 4. validate and plot results
for alignment in semi_repeat_free ahocorasick
do
	python3 ../scripts/compute_summary.py \
		-t 3 \
		--graph $inputgraph \
		--fastq output/sim_reads.fastq \
		--path output/sim_reads_path.nodes \
		--fasta output/sim_reads_path.fasta \
		--alignments output/best_${alignment}_alignments.gaf \
		--metrics output/metrics_${alignment}.mts &
done
wait $(jobs -p)

python3 ../scripts/compute_metrics.py \
	--output-name output/results \
	--summaries output/*.mts \
	--summaries-names $(basename -s ".mts" -a output/*.mts | sed 's/metrics_//')
