#!/usr/bin/env python

import argparse
import os
import sys
import io
import pandas as pd

from funpipe.gatk import gatk
import seaborn as sns
import matplotlib.pyplot as plt


stats = {
    'CountVariants': [
        'nVariantLoci', 'variantRatePerBp', 'nSNPs',
        'nInsertions', 'nDeletions', 'nHets', 'nHomRef', 'nHomVar',
        'hetHomRatio', 'insertionDeletionRatio'
    ],
    'TiTvVariantEvaluator': ['nTi', 'nTv', 'tiTvRatio'],
    'IndelSummary': ['SNP_to_indel_ratio']
}


def run_variant_eval(vcf, jar, fa, prefix, out_dir, RAM):
    prefix = os.path.join(out_dir, prefix)
    gatk_cmd = gatk(fa, jar, prefix)
    var_eval_tsv = gatk_cmd.variant_eval(vcf)
    return var_eval_tsv


def parse_variant_eval(eval):
    """ parse variantEval file
    :param eval: input eval file f
    :rtype
    """
    with open(eval, 'r') as fh:
        data = fh.read()
    tabs = data.rstrip().split('\n\n')

    meta_df = pd.DataFrame()
    for i in tabs:
        df = pd.read_csv(io.StringIO(i), comment='#', sep=r'\s+',
                         index_col='Sample')
        tab_name = df.columns[0]
        if tab_name in stats.keys():
            df = df[stats[tab_name]]
            meta_df = pd.concat([meta_df, df], axis=1)
    return meta_df


def parse_filter_geno_stat(file_geno_tsv):
    """ Load summary statistics from filterGenotypes.py """
    df = pd.read_csv(file_geno_tsv, sep='\t', header=0, index_col='Sample').T
    df.index.name = 'Sample'
    df.columns.name = None
    return df


def main(prefix, jar, out_dir, eval_tsv, filter_geno_stat, fa, RAM, vcf):
    if vcf:
        eval_tsv = run_variant_eval(vcf, jar, fa, prefix, out_dir, RAM)
        df = parse_variant_eval(eval_tsv)
    elif vcf:
        df = parse_variant_eval(eval_tsv)
    else:
        raise ValueError("Please input either an eval file or VCF file")

    if filter_geno_stat:
        filter_geno_df = parse_filter_geno_stat(filter_geno_stat)
        df = pd.concat([df, filter_geno_df], axis=1, join='inner')

    df.to_csv(os.path.join(out_dir, prefix+'.tsv.gz'), sep='\t',
              compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse output of GATK variant Eval')
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-p', '--prefix', help="Prefix of output file", required=True)
    required.add_argument('--jar', help='GATK jar')
    required.add_argument('--fa', help='reference fasta file')

    # optional arguments
    parser.add_argument(
        '-d', '--out_dir', default='.', help='Output Directory')
    parser.add_argument('-e', '--eval_tsv', help='Input eval file')
    parser.add_argument(
        '-f', '--filter_geno_tsv',
        help='filter summary statistics from filterGenotypes.py')
    parser.add_argument('--RAM', help='RAM', type=int, default=4)
    parser.add_argument('-v', '--vcf', help='Input vcf file')
    args = parser.parse_args()

    main(args.prefix, args.jar, args.out_dir, args.eval_tsv,
         args.filter_geno_tsv, args.fa, args.RAM, args.vcf)
