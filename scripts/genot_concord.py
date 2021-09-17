#!/usr/bin/env python3

import os
import sys
import argparse
from funpipe import gatk

fa = '/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta'
comp = ''
eval = '/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_cloud/output/CryptoD_test.genotype.filter.0608.vcf'


def main(eval, comp, ref):
    gatk_cmd = gatk(ref)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Compare concordance between two diploid VCFs'))
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-e', '--eval', required=True, help='VCF to be evaluated'
    )

    required.add_argument(
        '-c', '--comp', required=True, help='File as ground truth'
    )

    required.add_argument(
        '-r', '--ref', required=True, help='reference genome file'
    )

    args = parser.parse_args()
    main(args.eval, args.comp, args.ref)
