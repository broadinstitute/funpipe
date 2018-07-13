#!/usr/bin/env python3

from os.path import basename, splitext, join
import sys
import argparse
from pipeline import bam_depth, depth_per_window, sort_bam, cd, run, rm, clean_bam


def run_ploidy(out_dir, bam, faidx, no_sort):
    with cd(out_dir):
        base_prefix = splitext(basename(bam))[0]
        out_prefix = join(out_dir, base_prefix)
        bam_to_clean = bam
        if not no_sort:
            sorted_bam = sort_bam(bam, out_dir)
            bam_to_clean = sorted_bam
        cleaned_bam = clean_bam(bam_to_clean, out_prefix)
        pileup = bam_depth(cleaned_bam, out_prefix, idx=True)
        depth_per_window(pileup, base_prefix, faidx)
        # clean up working directory
        # if not no_sort:
        #     rm(sorted_bam)
        # rm(cleaned_bam)
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run ploidy analysis by analysing ')
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--bam', required=True, help='Input BAM path'
    )
    required.add_argument(
        '--faidx', required=True, help='Fasta index'
    )

    # optional arguments
    parser.add_argument(
        '-o', '--out_dir', help="Output file", default='.'
    )
    parser.add_argument(
        '-n', '--no_sort', help='Do not sort input bam', action='store_true',
        default=False
    )

    args = parser.parse_args()
    run_ploidy(args.out_dir, args.bam, args.faidx, args.no_sort)
