#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
import biolego.vcftools
import pandas as pd


comment_pattern = re.compile(r"^#")
qc_stats = pd.DataFrame(
    0, columns=samples,
    index=['filtered_'+x for x in ['GQ', 'AD', 'DP', 'Bi', 'all']])


def filter_genotypes(genotype_field, cutoffs):
    split_genotype = genotype_field.split(':')
    gq = (record.get_GQ(split_genotype[0],
                        index=sample_idx))
    init_gt = split_genotype[0]
    if ((split_genotype[0] in ['0', '0/0', '0|0'])
            and (cutoffs['keep_all_ref'])):
        (split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag) = (
            record.get_genotype(index=sample_idx, return_flags=True))
    elif (split_genotype[0] in ['0', '0/0', '0|0'] and
            gq == '0' and cutoffs['keep_GQ_0_refs']):
        (split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag) = (
            record.get_genotype(
                index=sample_idx,
                min_tot_dp=cutoffs['min_tot_DP'], return_flags=True))
    elif split_genotype[0] in ['0', '0/0', '0|0']:
        (split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag) = (
            record.get_genotype(
                index=sample_idx,
                min_gq=cutoffs['min_GQ'], min_tot_dp=cutoffs['min_tot_DP'],
                return_flags=True))
    else:
        (split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag) = (
            record.get_genotype(
                index=sample_idx,
                min_gq=cutoffs['min_GQ'], min_per_ad=cutoffs['min_AD'],
                min_tot_dp=cutoffs['min_tot_DP'],
                het_binom_p=cutoffs['het_binomial_p'],
                return_flags=True))
    if not (init_gt in ['.', './.']):
        if GQ_flag:
            qc_stats[sample]['filtered_GQ'] += 1
        if AD_flag:
            qc_stats[sample]['filtered_AD'] += 1
        if DP_flag:
            qc_stats[sample]['filtered_DP'] += 1
        if Bi_flag:
            qc_stats[sample]['filtered_Bi'] += 1
        if GQ_flag or AD_flag or DP_flag or Bi_flag:
            qc_stats[sample]['filtered_all'] += 1
    return str(":".join(split_genotype))


def get_site_info(record):
    site_info = [ str(record.get_chrom()),
                  str(record.get_pos()),
                  str(record.get_id()),
                  str(record.get_ref()),
                  str(record.get_alt_field()),
                  str(record.get_qual()),
                  str(record.get_filter()),
                  str(record.get_info()),
                  str(record.get_format())])
    return record


def main(input, outdir, prefix, cutoffs):

    header = vcfTools.VcfHeader(input)
    samples = header.get_samples()

    with open(infile) as vcf_file, open(prefix+'gt_ft.vcf', 'w') as gt_vcf:
        for vcf_line in vcf_file:
            if (re.search(comment_pattern, vcf_line)):
                gt_vcf.write(vcf_line)
            else:
                record = vcfTools.VcfRecord(vcf_line)
                site_info = get_site_info(record)
                for sample in samples:
                    sample_idx = header.get_sample_index(sample)
                    genotype_field = record.get_genotypes_field(record,
                                                                sample_idx)
                    gt_ft = filter_genotypes(genotype_fields, cutoffs)
                    site_info.append(gt_ft)
        gt_vcf.write("\t".join(site_info) + "\n")
    qc_stats.to_csv(prefix+'.stats.tsv', sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='filter vcf genotypes')
    parser.add_argument(
        '--het_binomial_p', default=False, type=float,
        help='filter hets that fail binomial test with given p-value')
    parser.add_argument(
        '--input', require=True, nargs='+',
        help='VCF filename')
    parser.add_argument(
        '--keep_all_ref', help='don\'t filter reference bases',
        action='store_true')
    parser.add_argument(
        '--keep_GQ_0_refs', action='store_true', # (keeps ref calls w RGQ/GQ=0)
        help='keep ref calls with GQ/RGQ 0 despite min_GQ')
    parser.add_argument(
        '--min_AD', type=float, default=0,
        help='minimum Allelic Depth')
    parser.add_argument(
        '--min_GQ', type=int, default=0,
        help='minimum GQ/RGQ to be kept')
    parser.add_argument(
        '--min_percent_alt_in_AD', type=float,
        help='min percent alt in AD to be kept for variants')
    parser.add_argument(
        '--min_total_DP', type=float, default=0,
        help='min total DP to be kept')
    parser.add_argument(
        '--outdir', help='output directory', default='./')
    parser.add_argument(
        '--prefix', help='output prefix', default="output")

    args = parser.parse_args()
    cutoffs = {
        "het_binomial_p": args.het_binomial_p,
        "keep_all_ref": args.keep_all_ref,
        "keep_GQ_0_refs": args.keep_GQ_0_refs,
        "min_AD": args.min_AD,
        "min_GQ": args.min_GQ,
        "min_percent_alt_in_AD":  args.min_percent_alt_in_AD,
        "min_total_DP": args.min_total_DP
    }
    main(args.input, args.outdir, args.prefix, cutoffs)
