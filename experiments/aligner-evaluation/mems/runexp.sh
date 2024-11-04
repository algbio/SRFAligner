#!/bin/bash
# remember to run `conda activate aligner-evaluation`
set -e
set -o pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

# executables' absolute paths/commands (make sure they work!)
BADREAD=$thisfolder/../../../tools/Badread/badread-runner.py
GRAPHCHAINER=$thisfolder/../../../tools/GraphChainer/bin/GraphChainer
graphaligner=$thisfolder/../../../tools/GraphAligner/bin/GraphAligner
srfaligner=$thisfolder/../../../SRFAligner
efgmemsaligner=$thisfolder/../../../efg-memsAligner
usrbintime=/usr/bin/time

# params
inputgraph=$thisfolder/../input/covid19_100_iEFG_simplified.gfa
coverage=1000
threads=1

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
	-m 0 \
	-i "N" \
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

$usrbintime $efgmemsaligner \
	-g $inputgraph \
	-f output/sim_reads.fastq \
	-a output/efg_mems_alignments.gaf \
	-w output \
	2>> output/runexp_log.txt >> output/runexp_log.txt

# 3. pick first alignment
for alignment in "semi_repeat_free_alignments.gaf" "graphaligner_alignments.gaf" "efg_mems_alignments.gaf"
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/$alignment > output/best_$alignment
done

# 4. validate and plot results
for alignment in semi_repeat_free graphaligner efg_mems
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
wait $(jobs -p)

python3 ../scripts/compute_metrics.py \
	--output-name output/results \
	--summaries output/*.mts \
	--summaries-names $(basename -s ".mts" -a output/*.mts | sed 's/metrics_//')
