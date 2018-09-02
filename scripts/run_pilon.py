#!/usr/bin/env python

import sys
import argparse
from funpipe.vcf import pilon, process_pilon_out


def main(inargs):
    args = parse_input_arg(inargs)
    with(cd(args.outdir)):
        fq1, fq2 = bam2fqs(args.bam, args.prefix, args.ram, args.picard_jar)
        realign_bam = bwa_align(args.fa, fq1, fq2, args.prefix)
        pilon(args.fa, realign_bam, args.prefix, args.ram, args.threads,
              args.pilon_jar)
        process_pilon_out(args.prefix+'.pilon.log', '.')


def parse_input_arg(args):
    parser = argparse.ArgumentParser()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '--prefix', help='output prefix', required=True)
    required.add_argument(
        '--bam', help='input bam', required=True)
    required.add_argument(
        '--fa', help='reference genome fasta file to evaluate', required=True)
    required.add_argument(
        '--outdir', help='output directory')

    # optional arguments
    parser.add_argument('--gff3', help='GFF3 annotation')
    parser.add_argument(
        '--ram', default=16, type=int, help='RAM usage of input file')
    parser.add_argument(
        '--threads', type=int, help='Number of threads for pilon', default=4
    )
    parser.add_argument(
        '--picard_jar', help = 'Picard jar',
        default='/seq/software/picard-public/current/picard.jar')
    parser.add_argument(
        '--snpeff_jar', help='Jar to snpeff',
        default='/cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.jar')
    parser.add_argument(
        '--pilon_jar', help='Pilon Jar',
        default='/gsap/assembly_analysis/apps/prod/pilon/lib/pilon-1.12.jar')
    parser.add_argument(
        '--snpeff_db', help='snpEff database'
    )
    return parser.parse_args()

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
