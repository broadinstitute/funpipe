#!/usr/bin/env python3

import os
import sys
import argparse
from pipeline import gatk

fa = '/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta'
comp = ''
eval = '/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_cloud/output/CryptoD_test.genotype.filter.0608.vcf'

gatk_cmd = gatk(fa)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(
#         description='')
#     # required arguments
#     required = parser.add_argument_group('required arguments')
#     required.add_argument(
#         '-i', '--input', required=True, help='Input file')
#
#     # optional arguments
#     parser.add_argument(
#         '-d', '--outdir', default='.', help='Output Directory'
#     )
#     parser.add_argument(
#         '-o', '--output', help="Output file",
#         default=sys.stdout, type=argparse.FileType('w'))
#
#     args = parser.parse_args()
