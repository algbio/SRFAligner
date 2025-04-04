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
srfaligner=$thisfolder/../../../SRFAligner
usrbintime=/usr/bin/time

# params
inputgraph=$thisfolder/../input/chr22_iEFG.gfa
coverage=30
threads=64

# 0. setup
mkdir output
echo -n > output/runexp_log.txt
set -a
set +a

# 1. simulate path and reads
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
	-f output/sim_reads.fastq \
	-a output/semi_repeat_free_alignments.gaf \
	2>> output/runexp_log.txt >> output/runexp_log.txt

$usrbintime $graphaligner \
	-t $threads \
	-x "vg" \
	-g $inputgraph \
	-f output/sim_reads.fastq \
	-a output/graphaligner_alignments.gaf \
	2>> output/runexp_log.txt >> output/runexp_log.txt

for o in 1 5 10 50 100
do
	$usrbintime $srfaligner \
		-t $threads \
		-g $inputgraph \
		-o $o \
		-f output/sim_reads.fastq \
		-a output/srf_edge_longest_${o}_alignments.gaf \
		2>> output/runexp_log.txt >> output/runexp_log.txt
done

# 3. pick first alignment
for alignment in "semi_repeat_free_alignments.gaf" "graphaligner_alignments.gaf"
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/$alignment > output/best_$alignment
done
for o in 1 5 10 50 100
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/srf_edge_longest_${o}_alignments.gaf > output/best_srf_edge_longest_${o}_alignments.gaf
done

# 4. validate and plot results
for alignment in semi_repeat_free graphaligner
do
	python3 ../scripts/compute_summary.py \
		-t 3 \
		--graph $inputgraph \
		--fastq output/sim_reads.fastq \
		--path output/sim_reads_path.nodes \
		--fasta output/sim_reads_path.fasta \
		--alignments output/best_${alignment}_alignments.gaf \
		--metrics output/metrics_${alignment}.mts
done
for o in 1 5 10 50 100
do
	python3 ../scripts/compute_summary.py \
		-t 3 \
		--graph $inputgraph \
		--fastq output/sim_reads.fastq \
		--path output/sim_reads_path.nodes \
		--fasta output/sim_reads_path.fasta \
		--alignments output/best_srf_edge_longest_${o}_alignments.gaf \
		--metrics output/metrics_srf_edge_longest_${o}.mts
done
wait $(jobs -p)

python3 ../scripts/compute_metrics.py \
	--output-name output/results \
	--summaries output/*.mts \
	--summaries-names $(basename -s ".mts" -a output/*.mts | sed 's/metrics_//')
