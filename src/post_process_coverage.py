#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from glob import glob
# import matplotlib.pyplot as plt


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


def get_contig_sets(all_contigs, subgenome):
    ''' Get all contig sets from each subgenome
    :param all_contigs: a set include all contigs
    :param subgenome: array of suffice in contigs for each subgenome
    '''
    contig_map = {}   # a hash map between contig map and their orders
    n_contig_set = [0] * len(subgenome)   # number of contigs in each subgenome
    for i in range(len(subgenome)):
        for j in all_contigs:
            if subgenome[i] in j:
                n_contig_set[i] += 1
                if i != 0:
                    n_sub_cont = (int(j.replace(subgenome[i], ''))
                                  + n_contig_set[i-1])
                    contig_map[j] = n_sub_cont
                else:
                    contig_map[j] = int(j.replace(subgenome[i], ''))
    return contig_map


def cov_fc_to_den(cov, contig_prefix, subgenome, split=False):
    ''' from coverage fold change file to generate density file needed for the
    downstream matlab ploting code
    :param cov: coverage data.frame
    :param contig_prefix: prefix of contig names
    :param subgenome: suffix of contig names, standing for different subgenomes
    :param split: split by subgenomes, not yet implemented
    :return: a datafram with coverage entries
    '''
    # remove unused columns
    den = cov.drop(['end0', 'id'], axis=1)
    den['chr'] = den.chr.str.replace(contig_prefix, '')
    if split:
        # subset to only a subgenome
        den = den[den.chr.str.contains(subgenome)]
    else:
        chrs = set(den.chr)
        contig_map = get_contig_sets(chrs, subgenome)

    # substitute contig names to numbers
    den.replace({'chr': contig_map}, inplace=True)
    den.sort_values(['chr', 'start0'], axis=0, inplace=True)
    den.drop(['start0'], axis=1, inplace=True)    # reformat to den file

    # add additional columns for den file
    den.insert(loc=1, column='id', value=range(1, den.shape[0]+1))
    den.insert(loc=2, column='dos1', value=1)
    den.insert(loc=3, column='dos2', value=1.5)
    den.insert(loc=4, column='dos3', value=2)
    den.insert(loc=5, column='dos4', value=3)
    return den


def output_den_tsv(prefix, den):
    ''' output density data structure to tsv file
    :param prefix: output Prefix
    :param den: density matrix
    '''
    # output den files
    with open(prefix+'.samples.tsv', 'w') as samples:
        for sample in den.columns[6:].tolist():
            samples.write(sample+'\n')
    den.to_csv(prefix+'.den', index=False, header=False, sep='\t')
    return 1


def post_process_coverage(fc_list, cov_tsv, prefix):
    """
    :param fc_list:  fold change list tsv
    :param cov_tsv: output file coverage tsv name
    :param prefix: output file prefix
    :param g_flags: subgenome flags
    """

    # prefix = 'batch1_75_AD'
    g_flags = ['_A', '_D']
    # combine coverage tsv
    cov = combine_cov_fc(input)
    # output merged table
    cov.to_csv(output, sep='\t', index=False)
    den = cov_fc_to_den(cov, 'chr', g_flags)
    output_den_tsv(prefix, den)
    # To do: add option to filter a specific subgenome
    # if args.g_flags is None:
    #     output_den_tsv(args.prefix, cov_fc_to_den(cov))
    # else:
    #     for flag in args.g_flags:
    #         # output sample IDs
    #         output_den_tsv(args.prefix+flag, cov_fc_to_den(cov, contain=flag))
    #
    return


def main(arguments):
    parser = argparse.ArgumentParser(
        description='PLOT COVERAGE ACROSS GENOME FOR MULTIPLE STRAINS')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input tsv lists'
    )
    required.add_argument(
        '-o', '--output', help="Output file",
        default=sys.stdout, type=argparse.FileType('w')
    )
    required.add_argument(
        '--faidx', help='fasta index file from samtools'
    )
    required.add_argument(
        '--depth', '-d', help='depth profile from samtools'
    )
    required.add_argument(
        '--prefix', '-p', help='output prefix'
    )
    # optional arguments
    parser.add_argument(
        '--g_flags', help='subgenome flag'
    )
    parser.add_argument('prefix', help='Prefix of output files')
    args = parser.parse_args(arguments)
    print(args)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
