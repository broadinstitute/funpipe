#!/usr/bin/env python
"""
Generate a genotype matrix using only loss of function variants from a list of
snpEff annotated VCF files.
"""
import sys
import re
import argparse
from funpipe import vcftools


comment_pattern = re.compile(r"^#")

stops = {}
amb_stops = {}
gene_names = []
annotations = {}
sample_list = []    # list of samples in all VCFs


def process_vcf(fof_line, ignore_indel_over, gof):
    fof_line = fof_line.rstrip()
    header = vcftools.VcfHeader(fof_line)
    caller = header.get_caller()
    samples = header.get_samples()
    per_file_sample_list = []
    for sample in samples:
        sample_list.append(sample)
        per_file_sample_list.append(sample)

    with open(fof_line, 'r') as vcf_file:
        for vcf_line in vcf_file:
            if not (re.search(comment_pattern, vcf_line)):
                record = vcftools.VcfRecord(vcf_line)
                is_passing = record.is_passing(caller)
                chrom = record.get_chrom()
                for sample in per_file_sample_list:
                    genotype = record.get_genotype(
                        index=header.get_sample_index(sample), min_gq=0)
                    variant_type = record.get_variant_type(caller, genotype)
                    alt = record.get_alt(genotype)
                    variant_annot = False
                    variant_impact = False
                    variant_feature = False
                    variant_effect = False
                    variant_length = 0
                    # parse snpeff results
                    if alt:
                        variant_annot = record.get_snpeff_annot(alt)
                        variant_impact = record.get_snpeff_impact(variant_annot)
                        variant_feature = record.get_snpeff_feature(variant_annot)
                        variant_effect = record.get_snpeff_effect(variant_annot)
                        variant_length = record.get_variant_length(genotype)

                    position = str(record.get_pos())

                    ## if syn/low and set to ignore syn/low
                    ## elsif passing snp or passing indel within size range
                    ## elsif ambig snp or amb indel within size range
                    is_stop = False
                    is_fs = False
                    is_gof = False
                    if variant_effect:
                        is_stop = re.search(r"stop_gained", variant_effect)
                        is_fs = re.search(r"frameshift_variant", variant_effect)
                        is_gof = re.search(r"stop_loss", variant_effect)

                    if (is_passing
                            and ((not gof and (is_stop or is_fs)) or (gof and is_gof))
                            and (variant_length <= ignore_indel_over)):
                        if variant_type in ['INSERTION', 'DELETION', 'SNP']:
                            if not variant_feature in stops:
                                stops[variant_feature] = {}
                                gene_names.append(variant_feature)
                            if not sample in stops[variant_feature]:
                                stops[variant_feature][sample] = 1
                            else:
                                stops[variant_feature][sample] += 1
                            if not variant_feature in annotations:
                                annotations[variant_feature] = [variant_annot]
                            elif not variant_annot in annotations[variant_feature]:
                                annotations[variant_feature].append(variant_annot)
                        elif variant_type == 'uncalled_ambiguous':
                            if not variant_feature in amb_stops:
                                amb_stops[variant_feature] = {}
                            if not sample in amb_stops[variant_feature]:
                                amb_stops[variant_feature][sample] = '-'


def output_lof_matrix():
    """ Output Loss-of-function genotype matrix """
    sorted_gene_names = sorted(gene_names)

    print("#feature", end="\t")
    print("\t".join(sorted_gene_names), end="")
    print("")

    print("#annotations", end="")
    for ind_gene in sorted_gene_names:
        print("\t", end="")
        try:
            print(",".join(annotations[ind_gene]), end="")
        except:
            print("NONE", end="")
    print("")

    print("reference", end="")
    for ind_gene in sorted_gene_names:
        print("\t0", end="")
    print("")

    for sample in sample_list:
        print(sample, end="")
        for ind_gene in sorted_gene_names:
            print("\t", end= "")
            try:
                if counts:
                    print(stops[ind_gene][sample], end="")
                elif stops[ind_gene][sample]:
                    print("1", end="")
            except:
                try:
                    if counts:
                        print(amb_stops[ind_gene][sample], end="")
                    elif amb_stops[ind_gene][sample]:
                        print("1", end="")
                except:
                    print("0", end="")
        print("")


def main(infile, ignore_indel_over, gof, counts):
    with open(infile, 'r') as fof:
        for fof_line in fof:
            sys.stderr.write("Searching " + fof_line)
            process_vcf(fof_line, ignore_indel_over, gof)
    output_lof_matrix()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', help='file of VCF filenames, must be annotated with SNPEff', type=str)
    parser.add_argument('--max_inframe_indel_size', help='maximum allowable in-frame indel size', type=int, default=12)
    parser.add_argument('--ignore_indel_over', help='ignore indels > <int> bp', type=int, default=100)
    parser.add_argument('--GOF', help='look for gain of function instead of of loss of function', action='store_true')
    parser.add_argument('--counts', help='print count of LOF mutations instead of binary 0/1', action='store_true')
    args = parser.parse_args()

    main(args.infile, args.ignore_indel_over, args.GOF, args.counts)
