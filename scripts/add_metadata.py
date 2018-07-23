#!/usr/bin/env python3

import os
import sys
import argparse


df = pd.read_csv('batch1_75_AD_SNPs_GQ50_AD08_DP10_stat.tsv', sep='\t',
                 index_col=0).T
meta = pd.read_csv('/cil/shed/sandboxes/xiaoli/fungal-pipeline/analysis/crypto/metadata/Summary_metadata_batch_5_6_1_Crypto_Sample_seroD.tsv',
                   sep='\t', index_col=1)
(pd.merge(df, meta, how='left', left_index='Sample', right_index='Sample')
 .to_csv('batch1_75_AD_samples_filter_stat.tsv', sep='\t', index=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Merge VCF QC stats with metadata')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input file')

    # optional arguments
    parser.add_argument(
        '-d', '--outdir', default='.', help='Output Directory'
    )
    parser.add_argument(
        '-o', '--output', help="Output file",
        default=sys.stdout, type=argparse.FileType('w'))

    args = parser.parse_args()
