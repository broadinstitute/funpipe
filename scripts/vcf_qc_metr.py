#!/usr/bin/env python

import argparse
import os
import sys
import io
import pandas as pd

stats = {
    'CountVariants': [
        'nVariantLoci', 'variantRatePerBp', 'nSNPs',
        'nInsertions', 'nDeletions', 'nHets', 'nHomRef', 'nHomVar',
        'hetHomRatio', 'insertionDeletionRatio'
    ],
    'TiTvVariantEvaluator': ['nTi', 'nTv', 'tiTvRatio'],
    'IndelSummary': ['SNP_to_indel_ratio']
}


def run_variant_eval(vcf):
    """ TO DO """
    return 1


def parse_variant_eval(eval, out_dir, out_tsv):
    """ parse variantEval file
    :param eval: input eval file f
    :param outdir:
    :param outfile: output file
    :rtype
    """
    with open(eval, 'r') as fh:
        data = fh.read()
    tabs = data.rstrip().split('\n\n')

    meta_df = pd.DataFrame()
    for i in tabs:
        df = pd.read_csv(io.StringIO(i), comment='#', sep='\s+',
                         index_col='Sample')
        tab_name = df.columns[0]
        if tab_name in stats.keys():
            df = df[stats[tab_name]]
            meta_df = pd.concat([meta_df, df], axis=1)
    meta_df.to_csv(os.path.join(out_dir, out_tsv), sep='\t',
                   compression='gzip')
    return meta_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse output of GATK variant Eval')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-e', '--eval', required=True, help='Input file')
    # optional arguments
    parser.add_argument(
        '-d', '--out_dir', default='.', help='Output Directory'
    )
    parser.add_argument(
        '-o', '--out_tsv', help="Output tsv file",
    )

    args = parser.parse_args()

    parse_variant_eval(args.eval, args.out_dir, args.out_tsv)
