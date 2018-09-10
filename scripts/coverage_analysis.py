#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from glob import glob
import matplotlib
import matplotlib.pyplot as plt
from distutils.version import LooseVersion
from funpipe.utils import run


# To do: add option to filter a specific subgenome
def combine_coverage_profiles(input):
    ''' Combine coverage fold changes profiles from samples
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
    print(' - Combine coverage profiles done.')
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


def coverage_to_density(cov, contig_prefix, subgenome, split=False):
    ''' From coverage fold change file to density file for the
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
    if LooseVersion(str(pd.__version__)) >= LooseVersion("0.17.0"):
        den.sort_values(['chr', 'start0'], axis=0, inplace=True)
    else:
        den.sort(['chr', 'start0'], axis=0, inplace=True)
    den.drop(['start0'], axis=1, inplace=True)    # reformat to den file
    den.insert(loc=1, column='id', value=range(1, den.shape[0]+1))
    return den


def output_density_tsv(prefix, den, legacy):
    ''' output density data structure to tsv file
    :param prefix: output Prefix
    :param den: density matrix
    :param legacy: if output tsv in legacy mode, compatible with matlab code
    '''
    out_den = prefix
    # output den files
    with open(prefix+'.samples.tsv', 'w') as samples:
        for sample in den.columns[6:].tolist():
            samples.write(sample+'\n')
    if legacy:
        out_den += '.den'
        # add additional columns for den file
        den.insert(loc=2, column='dos1', value=1)
        den.insert(loc=3, column='dos2', value=1.5)
        den.insert(loc=4, column='dos3', value=2)
        den.insert(loc=5, column='dos4', value=3)
        den.to_csv(out_den, index=False, header=False, sep='\t')
    else:
        out_den += '_cov_den.tsv'
        den.to_csv(out_den, index=False, sep='\t')
    return out_den


def cal_chr_percent(cov, min_cov):
    ''' Calculate proportion of reads coming from each chromosome
    :param cov: coverage data.frame
    :param min_cov: minimum coverage to include a window in coverage calculation
    :return data.frame: percentage of contig coverage.
    '''
    # filter coverage matrix to remove background coverage
    cov_ft = cov.iloc[:, 4:].applymap(lambda x: 0 if x < min_cov else x)
    cov_ft.insert(0, 'chr', cov.chr)
    # calculate percentage of coverage
    chr_pct_cov = {}
    chrs = sorted(list(set(cov.chr)))   # order of contigs to compute coverage
    for sample in cov_ft.columns.values[1:]:
        sample_cov = cov_ft[sample].sum()   # sample coverage
        chr_pct_cov[sample] = []
        for i in chrs:
            total_chr_cov = (cov_ft.loc[cov_ft.chr == i, sample].sum())
            chr_pct_cov[sample].append(total_chr_cov/sample_cov)
    chr_pct_cov_df = pd.DataFrame.from_dict(chr_pct_cov)
    chr_pct_cov_df.insert(0, 'contigs', chrs)
    return chr_pct_cov_df


def cal_subg_percent(cov, min_cov, subg):
    ''' Calculate proportion of reads from each subgenome
    :param cov: coverage dataframe
    :param min_cov: minimum coverage to include a window in the calculation
    :param subgneome: list that contains suffix of
    :return type: pandas dataframe
    '''
    cov_ft = cov.iloc[:, 4:].applymap(lambda x: 0 if x < min_cov else x)
    cov_ft.insert(0, 'chr', cov.chr)
    subg_pct_cov = {}
    for sample in cov_ft.columns.values[1:]:
        sample_cov = cov_ft[sample].sum()   # sample coverage
        subg_pct_cov[sample] = []
        for i in subg:
            subg_cov = (cov_ft.loc[cov_ft.chr.str.contains(i), sample].sum())
            subg_pct_cov[sample].append(subg_cov/sample_cov)
    subg_pct_cov_df = pd.DataFrame(subg_pct_cov)
    subg_pct_cov_df.insert(0, 'subg', subg)
    return subg_pct_cov_df


def coverage_plot(cov_tsv, prefix, color_csv, legacy):
    """ generate coverage plot
    :param cov_tsv: coverage profile list, first
    :param prefix: output prefix
    :param legacy: to do
    """
    cmd = ' '.join(['coverage_plot.R', '-f', cov_tsv, '-p', prefix,
                    '-c', color_csv])
    if legacy:
        cmd += ' -l'
    run(cmd)
    print(' - Finish generating coverage plot.')
    return 1


def main(cov_tsv, prefix, legacy, cutoff, no_plot, g_flags):
    cov_df = combine_coverage_profiles(cov_tsv)
    cov_df.to_csv(prefix+'.tsv', sep='\t', index=False)
    contig_pct_cov_df = cal_chr_percent(cov_df, cutoff)
    contig_pct_cov_df.to_csv(prefix+'.pct_cov.tsv', sep='\t', index=False)
    subg_pct_cov_df = cal_subg_percent(cov_df, )
    # calculate density
    den_df = coverage_to_density(cov_df, 'chr', g_flags)
    density_tsv = output_density_tsv(prefix, den_df, legacy)
    if not no_plot:
        coverage_plot(density_tsv, prefix, color_csv, legacy)
    print(' - Finish coverage analysis.')
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Aggregate coverage profiles of a sample set and generate'
        +' coverage plot per sample')
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--cov_tsv', required=True, help='A table-delimited file with '
        + 'first column sample ID, second column coverage profile path'
    )
    required.add_argument(
        '-p', '--prefix', help='output prefix', required=True
    )
#    required.add_argument(
#        '-c', '--color_csv', help='Color profile for each contig'
#    )

    # optional arguments
    parser.add_argument(
        '--no_plot', action='store_true', help='not generating coverage plots'
    )
    parser.add_argument(
        '--g_flags', help='subgenome flag', default=['_A', '_D']
    )
    parser.add_argument(
        '--min_cov', help=('minimum coverage per window for it to be included'
        ' in the analysis'), default=0
    )
    parser.add_argument(
        '--legacy', help='output density file in legacy mode, for matlab code',
        action='store_true'
    )
    args = parser.parse_args()

    main(
        args.cov_tsv, args.prefix, args.legacy, args.cutoff, args.no_plot,
        args.g_flags
    )
