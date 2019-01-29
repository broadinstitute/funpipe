#!/usr/bin/env python3

import os
import sys
import argparse
from biolego import phylo
from biolego import gatk
from biolego import analysis
from biolego.utils import cd


def phylo_analysis(input, prefix, max_amb, snp_only, filter_variants):
    analy = analysis(input, prefix, outdir, fasta=)
    if snp_only:
        gatks = gatk(analy, jar='')
        snp = gatks.select_snps()
    if filter_variants:
        filtered_vcf = gatks.filter_genotypes()
    phylos = phylo(analy)
    phylos.vcf_snp_to_fasta(max_amb)
    phylos.FastTreeDP()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform phylogenetic analysis with FreeTree')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input VCF', nargs='+')
    required.add_argument(
        '--fasta', required=True, help='input reference fasta file')
    # optional arguments
    parser.add_argument(
        '-s', '--snp_only', action='store_true',
        help='Whether to generate phylogenetics with only SNPs')
    parser.add_argument(
        '--max_amb', default=100000,
        help='maximum number of samples with ambiguous calls for a site to be'
             + ' included, recommended number of samples 10%')
    parser.add_argument('--ram', help='RAM usage for GATK')
    parser.add_argument(
        '--filter_variants', action='store_true',
        help='Whether filter genotypes before constructing phylogenetic tree')
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
        phylo_analysis(args.input, args.prefix, args.max_amb, args.snp_only)
