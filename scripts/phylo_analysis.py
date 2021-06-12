#!/usr/bin/env python3
""" Perform phylogenetic analysis with FreeTree """

import os
import sys
import argparse
from funpipe.phylo import vcf_snp_to_fasta, FastTreeDP
from funpipe.vcf import filter_variants
from funpipe.utils import cd


def phylo_analysis(input, prefix, max_amb):
    outvcf = filter_variants(input, prefix+'.vcf')
    fa = vcf_snp_to_fasta(outvcf, prefix, max_amb)
    FastTreeDP(fa, prefix)
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-v', '--vcf', required=True, help='Input vcf', nargs='+')
    required.add_argument(
        '-r', '--fasta', required=True, help='input reference fasta file')

    # optional arguments
    parser.add_argument(
        '--max_amb', default=10,
        help='maximum number of samples with ambiguous calls for a site to be'
             + ' included, recommended number of samples 10%')
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
