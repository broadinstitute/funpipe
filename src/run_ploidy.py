#!/usr/bin/env python3

import os
import sys
import argparse
from pipeline import bam_depth, depth_per_window, sort_bam, cd

def run_ploidy(bam, faidx):
    sorted_bam = sort_bam(bam)
    pileup = bam_depth(sorted_bam)
    depth_per_window(pileup, os.path.basename(bam), faidx)
    # rm(sorted_bam)

def main(arguments):
    parser = argparse.ArgumentParser(
        description='run ploidy analysis by analysing ')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--bam', required=True, help='Input file')
    required.add_argument(
        '--faidx', required=True, help='fasta index')
    # optional arguments
    parser.add_argument('--prefix', help='Prefix of output files')
    parser.add_argument(
        '-o', '--out_dir', help="Output file", default='.')

    args = parser.parse_args(arguments)
    cd(args.out_dir)
    run_ploidy(args.bam, args.faidx)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
