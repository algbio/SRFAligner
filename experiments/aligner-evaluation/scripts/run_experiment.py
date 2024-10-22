from os import getenv
from dotenv import load_dotenv
from subprocess import run
from argparse import ArgumentParser, RawTextHelpFormatter


def run_experiment(args):
    load_dotenv()

    srfaligner = getenv('SRFALIGNER')
    run(
        f'{srfaligner} -t {args.threads} -f {args.fastq} -g {args.graph} '
        f'-a {args.alignment}_srfaligner.gaf -w output'.split()
    )

#    graphaligner = getenv('GRAPHALIGNER')
#    run(
#        f'{graphaligner} -t {args.threads} -x vg -f {args.fastq} -g {args.graph} --verbose '
#        f'-a {args.alignment}_graphaligner.gaf'.split()
#    )

#    minigraph = getenv('MINIGRAPH')
#    run(
#        f'{minigraph} -t {args.threads} -c {args.graph} {args.fastq}'.split(),
#        stdout=open(f'{args.alignment}_minigraph.gaf', 'wb')
#    )

#    srfchainer = getenv('SRFCHAINER')
#    run(
#        f'{srfchainer} -t {args.threads} -f {args.fastq} -g {args.graph} '
#        f'-a {args.alignment}_srfchainer.gaf -w output'.split()
#    )

#    graphchainer = getenv('GRAPHCHAINER')
#    run(
#        f'{graphchainer} -t {args.threads} -f {args.fastq} -g {args.graph} '
#        f'-a {args.alignment}_graphchainer.gam '.split()
#    )

#    minichain = getenv('MINICHAIN')
#    run(
#        f'{minichain} -t {args.threads} -c {args.graph} {args.fastq}'.split(),
#        stdout=open(f'{args.alignment}_minichain.gaf', 'wb')
#    )


if __name__ == '__main__':

    parser = ArgumentParser(
        description='''
               Run aligners GraphAligner, minigraph, and srfaligner on the vg(?)/gfa and fastq files specified.            
            ''',
        formatter_class=RawTextHelpFormatter
    )

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-g', '--graph', type=str, help='Input vg/gfa file', required=True)
    requiredNamed.add_argument('-fq', '--fastq', type=str, help='Input fastq file', required=True)
    requiredNamed.add_argument(
        '-a', '--alignment', type=str, help='Output gam/gaf files (without extension)', required=True
    )

    parser.add_argument('-t', '--threads', type=int, help='Number of threads', default=30)

    run_experiment(parser.parse_args())
