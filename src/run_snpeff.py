#!/usr/bin/env python3

import os
import sys
import argparse
from pipeline import *

def main(arguments):
    parser = argparse.ArgumentParser(
        description='Run snpeff')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input_vcf', required=True, help='input vcf')
    required.add_argument(
        '-o', '--output_vcf', required=True, help='Output vcf')
    required.add_argument(
        '--genome_name', required=True,
        help='genome name in snpeff config file')
    required.add_argument(
        '-c', '--config', help='config file for snpeff', required=True)
    # optional arguments
    parser.add_argument('-m', '--ram', type=int, default=4, help='RAM usage')
    parser.add_argument('--snpeff_jar', help='jar file of snpeff',
        default='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db/snpEff.jar')

    args = parser.parse_args(arguments)

    snpeff(args.input_vcf, args.output_vcf, args.snpeff_jar, args.config,
        args.genome_name, args.ram)
    print(args)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
