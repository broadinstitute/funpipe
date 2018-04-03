#!/usr/bin/env python3

from os.path import basename, splitext, join
import sys
import argparse
from pipeline import bam_depth, depth_per_window, sort_bam, cd, run


def run_ploidy(out_dir, bam, faidx):
    with cd(args.out_dir):
        out_prefix = join(out_dir, splitext(basename(bam))[0])
        sorted_bam = sort_bam(bam)
        pileup = bam_depth(sorted_bam, out_prefix)
        depth_per_window(pileup, basename(bam), faidx)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run ploidy analysis by analysing ')
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--bam', required=True, help='Input BAM path')
    required.add_argument(
        '--faidx', required=True, help='Fasta index')

    # optional arguments
    parser.add_argument(
        '-o', '--out_dir', help="Output file", default='.')

    args = parser.parse_args()
    run_ploidy(args.out_dir, args.bam, args.faidx)
