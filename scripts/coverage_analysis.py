#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from glob import glob
import matplotlib
import matplotlib.pyplot as plt

# To do: add option to filter a specific subgenome
# if args.g_flags is None:
#     output_den_tsv(args.prefix, cov_fc_to_den(cov))
# else:
#     for flag in args.g_flags:
#         # output sample IDs
#         output_den_tsv(args.prefix+flag,
#                        cov_fc_to_den(cov, contain=flag))
#


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
    den.insert(loc=1, column='id', value=range(1, den.shape[0]+1))
    return den


def output_den_tsv(prefix, den, legacy):
    ''' output density data structure to tsv file
    :param prefix: output Prefix
    :param den: density matrix
    :param legacy: if output tsv in legacy mode, compatible with matlab code
    '''
    # output den files
    with open(prefix+'.samples.tsv', 'w') as samples:
        for sample in den.columns[6:].tolist():
            samples.write(sample+'\n')
    if legacy:
        # add additional columns for den file
        den.insert(loc=2, column='dos1', value=1)
        den.insert(loc=3, column='dos2', value=1.5)
        den.insert(loc=4, column='dos3', value=2)
        den.insert(loc=5, column='dos4', value=3)
        den.to_csv(prefix+'.den', index=False, header=False, sep='\t')
    else:
        den.to_csv(prefix+'_cov_den.tsv', index=False, sep='\t')
    return 1


def cal_subg_percent(cov, threhold=1):
    ''' calculate subgenome percentage for each chromosome
    :param cov: coverage data.frame
    :return data.frame: percentage of contig coverage.
    '''
    # subgenome map
    # subg_map = {
    #     'chr1_A': 'chr1_D',
    #     'chr2_A': 'chr2_D',
    #     'chr3_A': 'chr12_D',
    #     'chr4_A': 'chr13_D',
    #     'chr5_A': 'chr6_D',
    #     'chr6_A': 'chr7_D',
    #     'chr7_A': 'chr8_D',
    #     'chr8_A': 'chr9_D',
    #     'chr9_A': 'chr11_D',
    #     'chr10_A': 'chr12_D',
    #     'chr11_A': 'chr5_D',
    #     'chr12_A': 'chr15_D',
    #     'chr13_A': 'chr16_D',
    #     'chr14_A': 'chr10_D'
    # }
    # calculate proportion of reads coming from each chromosome
    # wether coverage above threshold
    cov_thres = cov.iloc[:, 4:].apply(lambda x: x >= 1)
    cov_thres.insert(0, 'chr', cov.contigs)
    chrs = sorted(list(set(cov.chr)))
    contig_pct_cov = {}
    for sample in cov_thres.columns.values[2:]:
        sample_cov = cov[sample].sum()
        contig_pct_cov[sample] = []
        for i in chrs:
            pass_thres = (cov_thres.loc[cov_thres.chr == i, sample].sum())
            contig_pct_cov[sample].append(pass_thres)
    contig_pct_cov_df = pd.DataFrame.from_dict(contig_pct_cov)
    contig_pct_cov_df.insert(0, 'contigs', chrs)
    return contig_pct_cov_df


def cmp_thres(val, threshold):
    if val >= threshold:
        return 1
    else:
        return 0


def coverage_analysis(fc_list, prefix, legacy):
    """
    :param fc_list:  fold change list tsv
    :param prefix: output file prefix
    :return
    """
    # prefix = 'batch1_75_AD'
    g_flags = ['_A', '_D']
    # combine coverage tsv
    cov = combine_cov_fc(fc_list)
    den = cov_fc_to_den(cov, 'chr', g_flags)
    contig_pct_cov_df = cal_subg_percent(cov)
    # output
    cov.to_csv(prefix+'.tsv', sep='\t', index=False)  # coverage table
    output_den_tsv(prefix, den, legacy)
    contig_pct_cov_df.to_csv(prefix+'.pct_cov.tsv', sep='\t', index=False)
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='PLOT COVERAGE ACROSS GENOME FOR MULTIPLE STRAINS')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--fc_list', required=True, help='Input tsv lists'
    )
    required.add_argument(
        '-p', '--prefix', help='output prefix', required=True
    )

    # optional arguments
    # parser.add_argument(
    #     '--g_flags', help='subgenome flag'
    # )
    parser.add_argument(
        '--legacy', help='output density file in legacy mode, for matlab code',
        action='store_true', default=False
    )
    args = parser.parse_args()

    print(args)
    coverage_analysis(args.fc_list, args.prefix, args.legacy)