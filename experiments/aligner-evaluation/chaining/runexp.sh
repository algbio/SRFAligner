#!/bin/bash
# remember to run `conda activate aligner-evaluation`
set -e
set -o pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

# executable's absolute paths/commands (make sure they work!)
BADREAD=$thisfolder/../../../tools/Badread/badread-runner.py
GRAPHCHAINER=$thisfolder/../../../tools/GraphChainer/bin/GraphChainer
srfchainer=$thisfolder/../../../SRFChainer
usrbintime=/usr/bin/time

# params
inputgraph=$thisfolder/../input/t2t-1KGP-phased-hapl-complete-250-heur-unsimplified.gfa
coverage=30
threads=64

# 0. setup
mkdir output
echo -n > output/runexp_log.txt
set -a
set +a

# 1. simulate path and reads
# comment the following 3 lines and uncomment the following ones to generate the reads (again)
ln -s $thisfolder/../semi-repeat-free/output/sim_reads.fastq output
ln -s $thisfolder/../semi-repeat-free/output/sim_reads_path.fasta output
ln -s $thisfolder/../semi-repeat-free/output/sim_reads_path.nodes output
#python3 ../scripts/generate_sim_reads.py \
#	--graph $inputgraph \
#	--fastq output/sim_reads.fastq \
#	--seed $(openssl enc -aes-256-ctr -pass pass:"semirepeatfree" -nosalt </dev/zero 2>/dev/null | shuf -i0-4294967295 -n 1 --random-source=/dev/stdin) \
#	--path output/sim_reads_path.nodes \
#	--fasta output/sim_reads_path.fasta \
#	--coverage $coverage \
#	2>> output/runexp_log.txt >> output/runexp_log.txt

# 2. run the aligners
$usrbintime $srfchainer \
	-t $threads \
	-g $inputgraph \
	-m 0 \
	-f output/sim_reads.fastq \
	-a output/semi_repeat_free_chain_alignments.gaf \
	2>> output/runexp_log.txt >> output/runexp_log.txt

for m in 1 2 3 4 5 6 7 8 9 10
do
	$usrbintime $srfchainer \
		-t $threads \
		-g $inputgraph \
		-m $m \
		-f output/sim_reads.fastq \
		-a output/srf_edge_${m}_chain_alignments.gaf \
		2>> output/runexp_log.txt >> output/runexp_log.txt
done

# 3. pick first alignment
for alignment in "semi_repeat_free_chain_alignments.gaf"
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/$alignment > output/best_$alignment
done
for m in 1 2 3 4 5 6 7 8 9 10
do
	awk '{if (found[$1] == "1") \
	          {} \
	      else
	          {found[$1]="1"; print}}' \
		output/srf_edge_${m}_chain_alignments.gaf > output/best_srf_edge_${m}_chain_alignments.gaf
done

# 4. validate and plot results
for alignment in semi_repeat_free_chain
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
for m in 1 2 3 4 5 6 7 8 9 10
do
	python3 ../scripts/compute_summary.py \
		-t 3 \
		--graph $inputgraph \
		--fastq output/sim_reads.fastq \
		--path output/sim_reads_path.nodes \
		--fasta output/sim_reads_path.fasta \
		--alignments output/best_srf_edge_${m}_chain_alignments.gaf \
		--metrics output/metrics_srf_edge_${m}_chain.mts &
done
wait $(jobs -p)

python3 ../scripts/compute_metrics.py \
	--output-name output/results \
	--summaries output/*.mts \
	--summaries-names $(basename -s ".mts" -a output/*.mts | sed 's/metrics_//')
