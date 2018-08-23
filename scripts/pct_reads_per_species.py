#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import os
import re


def reads_per_prot(infile):
    reads_per_prot = {}
    with open(infile, 'r') as input:
      for line in input:
        prot_id = line.strip().split('\t')[1].split('|')[1]
        if prot_id in reads_per_prot:
            reads_per_prot[prot_id] += 1
        else:
            reads_per_prot[prot_id] = 1
    return reads_per_prot


def prot_tax_map(taxonmap):
    prot_tax_map = {}
    with open(taxonmap, 'r') as map:
        for line in map:
            if '[' in line:
                record = re.split('[\[\]\|]', line.strip())
                prot_id, tax = record[1], record[4]
                prot_tax_map[prot_id] = tax
    return prot_tax_map


def reads_per_tax(prot_tax_map, reads_per_prot):
    reads_per_tax = {}
    total_reads = 0
    reads_per_tax['other'] = 0
    for prot_id, reads in reads_per_prot.items():
        total_reads += reads
        if prot_id in prot_tax_map:
            species = prot_tax_map[prot_id]
            if prot_id in reads_per_tax:
                reads_per_tax[species] += reads
            else:
                reads_per_tax[species] = reads
        else:
            reads_per_tax['other'] += reads
    return reads_per_tax, total_reads


def pct_reads_per_tax(reads_per_tax, total_reads):
    pct_reads_per_tax = {}
    for species, reads in reads_per_tax.items():
        pct = round(reads/total_reads * 100, 3)
        if pct >= 0.0001:
            pct_reads_per_tax[species] = pct
    return pct_reads_per_tax


def output(pct_reads_per_tax, prefix):
    (pd.DataFrame.from_dict(pct_reads_per_tax, orient='index')
        .sort_values(by=[0], ascending=False)
        .to_csv(prefix+'_pct_species.tsv', sep='\t', header=None))
    return 1


def main(taxonids, diamond_blastx_tsv, prefix):
    prot_tax_map_dict = prot_tax_map(taxonids)
    reads_per_prot_dict = reads_per_prot(diamond_blastx_tsv)
    (reads_per_tax_dict, total_reads) = reads_per_tax(prot_tax_map_dict, reads_per_prot_dict)
    pct_reads_per_tax_dict = pct_reads_per_tax(reads_per_tax_dict, total_reads)
    output(pct_reads_per_tax_dict, prefix)
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse diamond alignment results and derive percent of reads'
            + 'from each '
    )
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--diamond_blastx_tsv', required=True, help='Diamond output blast tabular format'
    )
    required.add_argument(
        '--taxonids', required=True,
        help='Path to mapping file that maps NCBI protein accession numbers to taxon ids'
    )
    # optional arguments
    parser.add_argument(
        '-d', '--out_dir', help='', default='.'
    )
    parser.add_argument(
        '-p', '--prefix', help='Output prefix', default='outfile'
    )

    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        raise ValueError(args.out_dir+"not exist.")

    main(args.taxonids, args.diamond_blastx_tsv,
         os.path.join(args.out_dir, args.prefix))
