#!/bin/bash
set -e
set -o pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

# executable's absolute paths/commands (make sure they work!)
bwa=$thisfolder/../../tools/bwa/bwa
efglocate=$thisfolder/../../tools/efg-locate/efg-locate
seqtk=$thisfolder/../../tools/seqtk/seqtk
usrbintime=/usr/bin/time

# params
inputgraph=$thisfolder/input/chr22_iEFG.gfa
inputreference=$thisfolder/input/chr22_uppercase.fasta
inputreads=$thisfolder/input/ERR1025645_sample05_1.fq.gz
threads=16

# 0. setup
mkdir output
echo -n > output/runexp_log.txt

echo "# 1. match in the chr22 iEFG (efg-locate)" >> output/runexp_log.txt
/usr/bin/time $efglocate \
	--reverse-complement \
	--threads $threads \
	$inputgraph \
	<(seqtk seq -A $inputreads) \
	output/efg_locate_matches.gaf \
	>> output/runexp_log.txt 2>> output/runexp_log.txt

echo "# 2. match in the chr22 reference (bwa)" >> output/runexp_log.txt
ln -s $inputreference $thisfolder/output/ref.fasta
/usr/bin/time $bwa index output/ref.fasta 2>> output/runexp_log.txt
/usr/bin/time $bwa aln -n 0 -k 0 -l 100 -t $threads \
	output/ref.fasta \
	$inputreads \
	> output/bwa_matches.sai 2>> output/runexp_log.txt
/usr/bin/time $bwa samse \
	output/ref.fasta \
	output/bwa_matches.sai \
	$inputreads \
	> output/bwa_matches.sam 2>> output/runexp_log.txt

echo "# 3. compute stats" >> output/runexp_log.txt
echo -n "efg-locate took" $(grep system output/runexp_log.txt | head -n 1 | cut -d' ' -f3 | cut -d'e' -f1) >> output/runexp_log.txt
echo " and matched" $(cut -f1 output/efg_locate_matches.gaf | uniq | sort | uniq | wc -l) "reads" >> output/runexp_log.txt

echo -n "bwa took" $(grep system output/runexp_log.txt | tail -n 3 | cut -d' ' -f3 | cut -d'e' -f1 | tr "\n" " ") >> output/runexp_log.txt
echo " and matched" $(cat output/bwa_matches.sam | awk '{if ($6 == "100M") {print}}' | cut -f1 | uniq | sort | uniq | wc -l) "reads" >> output/runexp_log.txt
