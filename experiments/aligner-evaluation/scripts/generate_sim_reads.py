from os import getenv
from dotenv import load_dotenv
from subprocess import run
from argparse import ArgumentParser, RawTextHelpFormatter


def generate_sim_reads(args):

    load_dotenv()
    badread = getenv('BADREAD')
    graphchainer = getenv('GRAPHCHAINER')

    run(
        f'{graphchainer} --generate-path --generate-path-seed {args.seed} -g {args.graph} -f {args.fasta} -a tmp.gam'
        .split()
    )

    run(f'mv tmp.gam {args.path}'.split())

    run(
        f'python {badread} simulate --identity 95,99,2.5 --seed {args.seed} --reference {args.fasta} '
        f'--quantity {args.coverage}x --length 15000,13000 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 '
        f'--chimeras 0'.split(),
        stdout=open(args.fastq, 'wb')
    )


if __name__ == '__main__':

    parser = ArgumentParser(
        description='''
            Generates simulated reads from a random path of an input vg/gfa file using the Badread simulator.
            Badread parameters are fixed to --identity 95,99,2.5 --length 15000,13000 --error_model nanopore2023
            --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0.
        ''',
        formatter_class=RawTextHelpFormatter
    )

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-g', '--graph', type=str, help='Input vg/gfa file', required=True)
    requiredNamed.add_argument('-fq', '--fastq', type=str, help='Output fastq file', required=True)

    parser.add_argument('-s', '--seed', type=int, help='Seed for random path generator and Badread', default=0)
    parser.add_argument(
        '-p', '--path', type=str, help='Output path file (node ids of selected path)', default='tmp.path'
    )
    parser.add_argument('-fa', '--fasta', type=str, help='Output fasta of original path', default='tmp.fasta')
    parser.add_argument('-c', '--coverage', type=int, help='Coverage value given to Badread', default=15)

    generate_sim_reads(parser.parse_args())

