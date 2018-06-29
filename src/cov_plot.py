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
        print(cov_list['Sample'][i])
        tab = (pd.read_csv(cov_list['Path'][i], sep='\t',
                           usecols=['chr', 'start0', 'end0', 'id', 'fc'])
               .rename(columns={'fc': cov_list['Sample'][i]}))
        if cov.empty:
            cov = tab
        else:
            cov = pd.merge(cov, tab, on=['chr', 'start0', 'end0', 'id'])
        print(cov.shape)
    return cov


def cov_fc_to_den(cov, contain=None):
    ''' from coverage fold change file to generate density file needed for the
    downstream matlab ploting code
    :param cov: coverage data.frame
    :param contain: characters defined subgenomes
    :return: a datafram with den entries
    '''
    # remove unused columns
    den = cov.drop(['end0', 'id'], axis=1)
    # subset to only a subgenome
    if contain is not None:
        den = den[den.chr.str.contains(contain)]
    # change chr to numbers and sort according to chromosomes
    den['chr'] = den.chr.str.replace(contain, '').str.replace('chr', '')
    den = den.sort(['chr', 'start0']).drop(['start0'], axis=1)
    # add additional columns for den file
    den.insert(loc=1, column='id', value=range(1, den.shape[0]+1))
    den.insert(loc=2, column='dos1', value=1)
    den.insert(loc=3, column='dos2', value=1.5)
    den.insert(loc=4, column='dos3', value=2)
    den.insert(loc=5, column='dos4', value=3)
    return den


def coverage_plot(cov):
    plt.show()
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

    prefix = 'batch1_75_AD'
    contain = ['_A', '_D']
    cov = combine_cov_fc(args.input)
    # output merged table
    cov.to_csv(args.output, sep='\t', index=False)
    for i in contain:
        # output sample IDs
        den = cov_fc_to_den(cov, contain=i)
        # samples list
        with open(prefix+i+'.samples.tsv', 'w') as samples:
            for sample in den.columns[6:].tolist():
                samples.write(sample+'\n')

        den.to_csv(prefix+i+'.den', index=False, header=False, sep='\t')
        # output den files

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
