#!/bin/bash
# sample some haplotypes from reference + VCF file, build MSA with vcf2multialign, and build iEFG with founderblockgraph
set -e
set -o pipefail

thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script

#
# executables
#
founderblockgraph=$thisfolder/../../tools/founderblockgraphs/founderblockgraph
vcf2multialign=vcf2multialign
bcftools=bcftools

#
# setup
#
threads=8
heuristicsubset=""

# parsing command line options
print_help()
{
   echo "usage: $0 [-f reference.fa] [-v variation.vcf] [-c chromosome] [-s nhapl] [-M Mrows] [-t threads]"
   echo "	-h --help:  show this screen"
   echo "	-f reference:  reference (FASTA format)                 used to generate MSA"
   echo "	-v variation:  variation (VCF format, possibly gzipped) used to generate MSA"
   echo "	-c chromosome: chromosome name                          used to generate MSA (see bcftools and vcf2multialign)"
   echo "	-s nhapl: sample the haplotype list and keep 'nhapl' random haplotypes"
   echo "	-M Mrows: set parameter --heuristic-subset='Mrows' for iEFG construction (see founderblockgraph)"
   echo "	-t threads: threads used in iEFG construction                            (see founderblockgraph)"
}

# https://stackoverflow.com/questions/12022592/how-can-i-use-long-options-with-the-bash-getopts-builtin
for arg in "$@"; do
	shift
	case "$arg" in
		'--help')                set -- "$@" '-h'   ;;
		*)                       set -- "$@" "$arg" ;;
	esac
done

argf=false ; argc=false ; argv=false ; args=false ; argM=false ;
OPTIND=1
while getopts "hf:v:c:r:s:M:t:" option; do
	case $option in
		f) # fasta + vcf input : fasta
			argf=true
			reference="$(realpath $OPTARG)" ;;
		v) # fasta + vcf input : vcf
			argv=true
			vcf="$(realpath $OPTARG)" ;;
		c) # fasta + vcf input : chromosome for bcftools/vcf2multialign
			argc=true
			chromosome="$OPTARG" ;;
		s) # number of samples
			args=true
			nhapl="$OPTARG" ;;
		M) # heuristic subset
			argM=true
			heuristicsubset="--heuristic-subset $OPTARG" ;;
		t) # threads
			argt=true
			threads="$OPTARG" ;;
		h) # display help
			print_help
			exit;;
		\?) # invalid option
			echo "Error: Invalid option"
			exit;;
	esac
done
shift $(expr $OPTIND - 1) # remove options from positional parameters

if [ "$argf" = false ] || [ "$argc" = false ] || [ "$argv" = false ] || [ "$args" = false ]
then
	print_help
	exit
fi

outputfolder=$thisfolder/output
mkdir $outputfolder
if [[ $? -gt 0 ]] ; then echo "Output directory $outputfolder already exists!" ; exit 1 ; fi
log=$outputfolder/log.txt
stats=$outputfolder/stats.txt
cd $outputfolder

#
# randomness source
#
# https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

#
# 1. sampling the vcf
#
echo -n "Sampling the vcf..."
# TODO I am expecting the haplotypes to be after the 9th field specified by the VCF header, is this always correct?
$bcftools query -l $vcf > haplotypes
cat haplotypes | shuf -n $nhapl --random-source=<(get_seeded_random "semi-repeat-free") > sampled_haplotypes
$bcftools view --regions $chromosome --samples-file sampled_haplotypes $vcf > sampled_haplotypes.vcf
# FIX for the T2T 1KGP data and vcf2multialign, see https://github.com/tsnorri/vcf2multialign/issues/5
sed -i 's/e+06//g' sampled_haplotypes.vcf
sed -i 's/INFO=<ID=AC,Number=1/INFO=<ID=AC,Number=A/g' sampled_haplotypes.vcf
echo " done."

#
# 2. vcf -> MSA
#
echo -n "Computing the MSA..."
/usr/bin/time $vcf2multialign \
	--founder-sequences=50 \
	--input-reference=$reference \
	--input-variants=sampled_haplotypes.vcf \
	--chromosome $chromosome \
	--output-graph variant.graph >> $log 2>> $log

/usr/bin/time $vcf2multialign \
	--input-reference=$reference \
	--input-graph variant.graph \
	-H \
	-s sampled_haplotypes.a2m >> $log 2>> $log
echo " done."

#
# 3. MSA -> iEFG
#
echo -n "Building the indexable Elastic Founder Graph..."
/usr/bin/time $founderblockgraph --ignore-chars="N" --output-paths --threads=$threads --input=sampled_haplotypes.a2m --output=efg-unsimplified.gfa $heuristicsubset >> $log 2>> $log
echo " done."
