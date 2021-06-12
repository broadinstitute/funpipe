#!/usr/bin/env python3
import argparse
import os
from funpipe import bam.fastqc, picard


def bam_fast_qc(bam, ref_fa, out_dir, prefix):
    ''' use FastQc for a BAM '''
    picard_cmd = picard()
    (fq1, fq2) = picard_cmd.bam2fqs(bam, os.path.join(out_dir, prefix))
    fastqc(bam, fq1, fq2, out_dir)
    print(' - bam_fast_qc done.')
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Exam quality of BAMs'
    )
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--bam', required=True, help='Input BAM file'
    )
    required.add_argument(
        '-f', '--ref_fa', required=True, help='Reference genome fasta file'
    )
    # optional arguments
    parser.add_argument(
        '-d', '--out_dir', help='output directory', default='.'
    )
    parser.add_argument(
        '-p', '--prefix', help='Output prefix', default='outfile'
    )
    parser.add_argument(
        '-j', '--picard_jar', help='path to picard_jar',
        default='/seq/software/picard/1.853/bin/picard.jar'
    )

    args = parser.parse_args()
    bam_fast_qc(args.bam, args.ref_fa, args.out_dir, args.prefix)
