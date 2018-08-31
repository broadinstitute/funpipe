#!/usr/bin/env python3

import os
import sys
import argparse
from funpipe import cd, vcf_snp_to_fasta, FastTreeDP, filter_variants


def phylo_analysis(input, prefix, max_amb):
    outvcf = filter_variants(input, prefix+'.vcf')
    fa = vcf_snp_to_fasta(outvcf, prefix, max_amb)
    FastTreeDP(fa, prefix)
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform phylogenetic analysis with FreeTree')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input vcf', nargs='+')
    required.add_argument(
        '--max_amb', required=True,
        help='maximum number of samples with ambiguous calls for a site to be'
             + ' included, recommended number of samples 10%')
    required.add_argument(
        '--fasta', required=True, help='input reference fasta file')
    # optional arguments
    parser.add_argument('--ram', help='RAM usage for GATK')
    parser.add_argument(
        '-d', '--outdir', default='.', help='Output Directory')
    parser.add_argument(
        '-p', '--prefix', default='Output', help='Prefix of output file')
    args = parser.parse_args()

    with cd(args.outdir):
        # if args.input > 1:
        #     if fa is None:
        #         raise ValueError('Please specify fasta file.\n')
        #     else:
        #         (gatk(fa, RAM=args.ram)
        #          .combineVar())
        phylo_analysis(args.input, args.prefix, args.max_amb)
