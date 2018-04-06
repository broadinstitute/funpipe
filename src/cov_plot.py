#!/usr/bin/env python3

import os
import sys
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob


def combine_cov_fc(input):
    ''' combine coverage fold changes profiles from samples
    :param input: input file containing sample and
    :return Pandas DataFrame
    '''
    cov_list = pd.read_csv(
        input, sep='\t', header=None, names=('Sample', 'Path'))
    cov = pd.DataFrame()
    for i in range(0, cov_list.shape[0]):
        tab = (pd.read_csv(cov_list['Path'][i], sep='\t',
                           usecols=['chr', 'start0', 'end0', 'id', 'fc'])
               .rename(columns={'fc': cov_list['Sample'][i]}))
        if cov.empty:
            cov = tab
        else:
            cov = pd.merge(cov, tab, on=['chr', 'start0', 'end0', 'id'])
    return cov


def coverage_plot(cov):

    return


def main(arguments):
    parser = argparse.ArgumentParser(
        description='PLOT COVERAGE ACROSS GENOME FOR MULTIPLE STRAINS')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input tsv lists')
    required.add_argument(
        '-o', '--output', help="Output file",
        default=sys.stdout, type=argparse.FileType('w'))
    required.add_argument(
        '--faidx', help='fasta index file from samtools')
    required.add_argument(
        '--depth', '-d', help='depth profile from samtools')
    # optional arguments
    parser.add_argument('prefix', help='Prefix of output files')

    args = parser.parse_args(arguments)
    print(args)

    cov = combine_cov(input)



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
