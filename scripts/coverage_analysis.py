#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
from distutils.version import LooseVersion
from funpipe.utils import run


def combine_coverage_profiles(input, prefix=None):
    ''' Combine coverage fold changes profiles from samples
    :param input: input file containing sample and
    :param perfix: optional, when provided, will output cov_df to a tsv file.
    :return Pandas DataFrame
    '''
    cov_list = pd.read_csv(
        input, sep='\t', header=None, names=('Sample', 'Path'))
    cov_df = pd.DataFrame()
    for i in range(0, cov_list.shape[0]):
        tab = (pd.read_csv(cov_list['Path'][i], sep='\t',
                           usecols=['chr', 'start0', 'end0', 'id', 'fc'])
               .rename(columns={'fc': cov_list['Sample'][i]}))
        if cov.empty:
            cov_df = tab
        else:
            cov_df = pd.merge(cov_df, tab, on=['chr', 'start0', 'end0', 'id'])
    if prefix is not None:
        cov_df.to_csv(prefix+'_cov.tsv', sep='\t', index=False)
    print(' - Combine coverage profiles done.')
    return cov_df


def get_contig_sets(all_contigs, subgenome):
    ''' Get all contig sets from each subgenome
    :param all_contigs: a set include all contigs
    :param subgenome: array of suffice in contigs for each subgenome
    '''
    contig_map = {}  # a hash map between contig map and their orders
    n_contig_set = [0] * len(subgenome)  # number of contigs in each subgenome
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


def coverage_to_density(cov, contig_prefix, subgenome, split=False,
                        prefix=None, legacy=False):
    ''' From coverage fold change file to density file for the
    downstream matlab ploting code
    :param cov: coverage data.frame
    :param contig_prefix: prefix of contig names
    :param subgenome: suffix of contig names, standing for different subgenomes
    :param split: whether to produce a density file from a subgenome
    :param prefix: optional, when provided, will output cov_df to a tsv file.
    :param legacy: whether output tsv in legacy mode, compatible with matlab
                   code
    :return: a datafram with coverage entries
    '''
    # remove unused columns
    den = cov.drop(['end0', 'id'], axis=1)
    den['chr'] = den.chr.str.replace(contig_prefix, '')
    if split:  # subset to only a subgenome
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
    if prefix is not None:
        density_tsv = output_density_tsv(prefix, den, legacy)
    return den


def output_density_tsv(prefix, den, legacy):
    ''' output density data structure to tsv file
    :param prefix: output Prefix
    :param den: density matrix
    :param legacy: output tsv in legacy mode, compatible with matlab code
    :rtype pd dataframe
    '''
    out_den = prefix
    # output den files
    # with open(prefix+'.samples.tsv', 'w') as samples:
    #     for sample in den.columns[6:].tolist():
    #         samples.write(sample+'\n')
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


def filter_cov(cov, min_cov):
    """ Filter coverage dataframe
    :param cov: coverage data.frame
    :param min_cov: minimum coverage to include a window in analysis
    """
    cov_ft = cov.iloc[:, 4:].applymap(lambda x: 0 if x < min_cov else x)
    cov_ft.insert(0, 'chr', cov.chr)
    return cov_ft


def cal_chr_percent(cov, prefix=None):
    ''' Calculate proportion of reads coming from each chromosome
    :param cov: coverage data.frame
    :param prefix: optional, when provided, will output cov_df to a tsv file.
    :rtype data.frame: percentage of contig coverage.
    '''
    # filter coverage matrix to remove background coverage
    # calculate percentage of coverage
    chr_pct_cov = {}
    chrs = sorted(list(set(cov.chr)))  # order contigs for pct
    for sample in cov.columns.values[1:]:
        sample_cov = cov[sample].sum()   # sample coverage
        chr_pct_cov[sample] = []
        for i in chrs:
            total_chr_cov = (cov.loc[cov.chr == i, sample].sum())
            chr_pct_cov[sample].append(round(total_chr_cov/sample_cov, 4))
    chr_pct_cov_df = pd.DataFrame.from_dict(chr_pct_cov)
    chr_pct_cov_df.insert(0, 'contigs', chrs)
    if prefix is not None:
        chr_pct_cov_df.to_csv(prefix+'_pct_cov.tsv', sep='\t', index=False)
    print(" - Reads from each chromosome is calcuated.")
    return chr_pct_cov_df


def chr_coverage(cov, prefix=None):
    """ Calculate percentage of chromosomes covered by reads, by calculating
    proportion of bins passing minimum coverage threshold
    :param cov: coverage data.frame
    :param prefix: optional, when provided, will output cov_df to a tsv file.
    :rtype: pandas dataframe
    """
    cov_bin = cov.iloc[:, 1:].applymap(lambda x: 1 if x > 0 else 0)
    cov_bin.insert(0, 'chr', cov.chr)
    chrs = sorted(list(set(cov.chr)))   # order of contigs to compute coverage
    chr_cov = {}
    for sample in cov_bin.columns.values[1:]:
        chr_cov[sample] = []
        for i in chrs:
            n_total_bin = cov_bin.loc[cov_bin.chr == i, sample].shape[0]
            n_covered_bin = cov_bin.loc[cov_bin.chr == i, sample].sum()
            chr_cov[sample].append(round(n_covered_bin/n_total_bin, 4))
    chr_cov_df = pd.DataFrame.from_dict(chr_cov)
    chr_cov_df.insert(0, 'chr', chrs)
    if prefix is not None:
        chr_cov_df.to_csv(prefix+'_chr_pct_cov.tsv', sep='\t', index=False)
    return chr_cov_df


def cal_subg_percent(cov, subg, prefix=None):
    ''' Calculate proportion of reads from each subgenome
    :param cov: coverage dataframe
    :param subgneome: list that contains suffix of
    :param prefix: optional, when provided, will output cov_df to a tsv file.
    :return type: pandas dataframe
    '''
    subg_pct_cov = {}
    for sample in cov.columns.values[1:]:
        sample_cov = cov[sample].sum()   # sample coverage
        subg_pct_cov[sample] = []
        for i in subg:
            subg_cov = (cov.loc[cov.chr.str.contains(i), sample].sum())
            subg_pct_cov[sample].append(round(subg_cov/sample_cov, 4))
    subg_pct_cov_df = pd.DataFrame(subg_pct_cov)
    subg_pct_cov_df.insert(0, 'subg', subg)
    if prefix is not None:
        subg_pct_cov_df.to_csv(prefix+'_subg_pct.tsv', sep='\t', index=False)
    return subg_pct_cov_df


def subg_barplot(df, prefix):
    """ Plot subgenome proportion
    :param df: data frame
    :param prefix: output prefix
    """
    sns.set(style="whitegrid")
    # Initialize the matplotlib figure
    f, ax = plt.subplots(figsize=(5, df.shape[0]*0.2))

    sns.set_color_codes("pastel")
    subp_plt = sns.barplot(x='_D', y='samples', data=df,
                           label="Pct_D", color="g")
    subp_fig = subp_plt.get_figure()
    subp_fig.subplots_adjust(left=0.5)
    subp_fig.savefig(prefix+'_subg_pct.png')
    ax.set(xlim=(0, 1))
    sns.despine(left=True, bottom=True)


def coverage_scatterplot(cov_tsv, prefix, color_csv, legacy):
    """ generate coverage plot
    :param cov_tsv: coverage profile list, first
    :param prefix: output prefix
    :param legacy: if output tsv in legacy mode, compatible with matlab code
    """
    cmd = ' '.join(['coverage_scatterplot.R', '-i', cov_tsv, '-p', prefix,
                    '-c', color_csv])
    if legacy:
        cmd += ' -l'
    run(cmd)
    print(' - Finish generating coverage plot.')
    return 1


def cal_frac_aneu(ploidy, ploidy_list):
    """
    Calculate percentage of aneiploidy for each ploidy list

    Examples
    --------

    >>> ploidy = [0, 1, 2, 4, 4]

    >>> ploidy_list = [0, 1, 2, 3, 4]

    >>> cal_frac_aneu(ploidy, ploidy_list)

    [0.2, 0.2, 0.2, 0, 0.4]

    Parameters
    ----------

    ploidy :
        a list of ploidy

    ploidy_list :
        a list of ploidy

    Returns
    -------

    a list of ploidy fractions

    """
    total = len(ploidy)
    counts = collections.Counter(ploidy)
    frac = []
    for dos in ploidy_list:
        if counts[dos]:
            frac.append(round(counts[dos]/total, 2))
        else:
            frac.append(0)
    return frac


def pct_aneuploidy(cov_df, max_ploidy=4, prefix=None):
    """ Calculate percentage of aneuploidy per chromosome for each sample

    Notes
    -----
    Ploidy of each sliding window was rounded, and ploidy greater than
    max_ploidy was set to max_ploidy.

    Examples
    --------

    >>> cov_df = pd.DataFrame({
            'chr': ['chr1_A', 'chr1_A', 'chr2_A', 'chr2_A'],
            'sample0': [0.1, 1.2, 2.6, 3],
            'sample1': [3.9, 4.6, 2.1, 1.4]
        }, columns=['chr', 'sample0', 'sample1'])

    >>> pct_aneuploidy(cov_df)

        sample  ploidy  chr1_A  chr2_A
    0  sample0       0     0.5     0.0
    1  sample0       1     0.5     0.0
    2  sample0       2     0.0     0.0
    3  sample0       3     0.0     1.0
    4  sample0       4     0.0     0.0
    5  sample1       0     0.0     0.0
    6  sample1       1     0.0     0.5
    7  sample1       2     0.0     0.5
    8  sample1       3     0.0     0.0
    9  sample1       4     1.0     0.0

    Parameters
    ----------
    cov_df : :obj:`dataframe`
        Input coverage dataframe

    max_ploidy : `int` maximum ploidy in this
    Returns
    -------
    aneu_df : :obj:`dataframe`
        Percent of aneuploidy for each chromosome

    """
    chrs = cov_df['chr'].unique().tolist()
    samples = cov_df.columns[1:].tolist()
    ploidy_list = list(range(0, max_ploidy+1))
    # round coverage per window capped at max_ploidy
    round_df = (cov_df.drop(['chr'], axis=1).round()
                .apply(
                    lambda x: [y if y < max_ploidy else max_ploidy for y in x])
                )
    round_df.insert(0, 'chr', cov_df['chr'])
    # generate aneuploidy data frame
    aneu_df = pd.DataFrame({
        'sample': [val for val in samples for _ in range(0, len(ploidy_list))],
        'ploidy': ploidy_list * len(samples)
    })
    for chr in chrs:
        sample_frac = []
        for sample in samples:
            ploidy = round_df[round_df.chr == chr][sample].tolist()
            print(ploidy)
            sample_frac += cal_frac_aneu(ploidy, ploidy_list)
            print(sample_frac)
        chr_frac_df = pd.DataFrame({chr: sample_frac})
        aneu_df = pd.concat([aneu_df, chr_frac_df], axis=1, sort=False)

    if prefix is not None:
        aneu_df.to_csv(prefix+'_aneu_frac.tsv', sep='\t', index=False)
    return aneu_df


def chr_cov_plot(cov_df, prefix):
    """ Pct chr coverage plot """
    df.set_index('chr', inplace=True)
    f, ax = plt.subplots(figsize=(df.shape[1]*0.5, df.shape[0]*0.5))
    if cluster:
        print(" - Cluster coverage profile.")
        chr_pct_fig = sns.clustermap(df, row_cluster=True, col_cluster=True)
    else:
        print(" - Produce heatmap from coverage profile.")
        chr_pct_fig = sns.heatmap(df, vmin=0, vmax=1,
                                  cmap="YlGnBu").get_figure()
    chr_pct_fig.savefig(prefix+'_chr_pct_cov.png')
    return 1


def main(cov_tsv, prefix, legacy, min_cov, no_plot, g_flags, color_csv,
         no_cluster):
    print('Start to process '+prefix)
    cov_df = combine_coverage_profiles(cov_tsv, prefix)
    cov_ft_df = filter_cov(cov_df, min_cov)
    pct_cov_df = cal_chr_percent(cov_ft_df, prefix)
    chr_pct_cov_df = chr_coverage(cov_ft_df, prefix)
    aneu_frac = pct_aneuploidy(cov_ft_df, prefix=prefix)
    chr_cov_heatmap(chr_pct_cov_df, prefix)
    subg_pct_cov_df = cal_subg_percent(cov_ft_df, g_flags, prefix)

    subg_df_t = (subg_pct_cov_df.set_index('subg').transpose()
                 .rename_axis('samples').rename_axis(None, 1).reset_index())
    subg_barplot(subg_df_t, prefix)
    # calculate density
    density_tsv = coverage_to_density(cov_df, 'chr', g_flags, prefix, legacy,
                                      split=False)
    coverage_scatterplot(density_tsv, prefix, color_csv, legacy)
    print(' - Done.')
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Aggregate coverage profiles of a sample set and generate'
                     ' coverage plot per sample')
    )
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--cov_tsv', required=True, help='A table-delimited file with '
        + 'first column sample ID, second column coverage profile path'
    )
    required.add_argument(
        '-p', '--prefix', help='output prefix', required=True
    )
    required.add_argument(
       '-c', '--color_csv', help='Color profile for each contig', required=True
    )

    # optional arguments
    parser.add_argument(
        '--g_flags', help='subgenome flag', default=['_A', '_D']
    )
    parser.add_argument(
        '--min_cov', help=('minimum coverage per window for it to be included'
                           ' in the analysis'), default=0, type=float
    )
    parser.add_argument(
        '--legacy', help='output density file in legacy mode, for matlab code',
        action='store_true'
    )
    parser.add_argument(
        '--no_cluster', action='store_true',
        help='Perform hierachical clustering for the pct chr coverage plot'
    )
    parser.add_argument(
        '--split', action='store_true',
        help='whether present only a subgenome '
    )
    args = parser.parse_args()

    main(
        args.cov_tsv, args.prefix, args.legacy, args.min_cov, args.no_plot,
        args.g_flags, args.color_csv, args.no_cluster
    )
