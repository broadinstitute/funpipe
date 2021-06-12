#!/usr/bin/env python3
""" Produce fasta files corresponding from SNP only VCFs

Note that this is a piece of legacy code from the Broad's fungal group, and is
provided as is, without optimizing the performance, code style and
documentation etc. Only minimal changes were introduced to ensure compatibility
with the funpipe package.

"""
import sys
import re
import argparse
from funpipe import vcftools
import gzip


parser = argparse.ArgumentParser()
parser.add_argument(
    'infile',
    help='a txt file where each line is a path to a vcf', type=str)
parser.add_argument(
    '--max_amb_samples',
    help=('maximum number of samples with ambiguous calls for a site to be '
          'included'), type=int, default=10)
args = parser.parse_args()

infile = args.infile
max_amb = args.max_amb_samples  # 100000 as default in legacy code


# Initialize variables
ref_bases = {}
alt_bases = {}
passed_snp_positions = {}
sample_list = []
amb_pos_counts = {}
sample_translations = {}

line_number = 1
site_in_span_dels = 0

# process each vcf file
fof = open(infile)
for vcf_name in fof:
    if vcf_name[-2:] == 'gz':
        raise ValueError((
            "This legacy code does not support compressed VCFs. Please "
            "decompress all your bgzipped VCFs and re-ran this analysis."))

    sys.stderr.write("Searching " + vcf_name)
    # process individual VCF
    vcf_name = vcf_name.rstrip()
    header = vcftools.VcfHeader(vcf_name)
    caller = header.get_caller()
    samples = header.get_samples()
    contigs = header.get_contigs()

    if line_number == 1:
        for contig in contigs:
            passed_snp_positions[contig] = {}
            amb_pos_counts[contig] = {}

    if samples == ['SAMPLE']:
        samples = [vcf_name]
        sample_translations[vcf_name] = 'SAMPLE'
        sys.stderr.write("No sample name in " + vcf_name +
                         ", using file name.\n")

    per_file_sample_list = []
    for sample in samples:
        sample_list.append(sample)
        per_file_sample_list.append(sample)

    with open(vcf_name) as vcffh:
        for vcf_line in vcffh:
            if vcf_line[0] != '#':
                record = vcftools.VcfRecord(vcf_line)
                pass_or_fail = record.is_passing(caller)
                # skip sites with '*' marks in alt field
                if '*' in record.get_alt_field():
                    site_in_span_dels += 1
                    continue
                for sample in per_file_sample_list:
                    translation = sample
                    try:
                        translation = sample_translations[sample]
                    except:
                        pass
                    genotype = (record.get_genotype(
                        index=header.get_sample_index(translation), min_gq=0))
                    variant_type = record.get_variant_type(caller, genotype)
                    if pass_or_fail and not variant_type:
                        pass
                    elif pass_or_fail and variant_type == 'SNP':
                        chr = record.get_chrom()
                        pos = int(record.get_pos())

                        if sample not in alt_bases:
                            alt_bases[sample] = {}
                        if chr not in alt_bases[sample]:
                            alt_bases[sample][chr] = {}

                        alt_bases[sample][chr][pos] = record.get_alt(genotype)
                        if chr not in ref_bases:
                            ref_bases[chr] = {}
                        ref_bases[chr][pos] = record.get_ref()

                        if alt_bases[sample][chr][pos] != 'N':
                            passed_snp_positions[chr][pos] = True
                    else:
                        chr = record.get_chrom()
                        pos = int(record.get_pos())
                        if sample not in alt_bases:
                            alt_bases[sample] = {}
                        if chr not in alt_bases[sample]:
                            alt_bases[sample][chr] = {}
                        alt_bases[sample][chr][pos] = 'N'

                        if pos not in amb_pos_counts[chr]:
                            amb_pos_counts[chr][pos] = 1
                        else:
                            amb_pos_counts[chr][pos] += 1
    line_number += 1
fof.close()


sys.stderr.write(
    str(site_in_span_dels)
    + " sites filtered due to locating in spanning deletions.\n")


# write fasta files
sorted_chrs = sorted(passed_snp_positions.keys())

print(">reference")
sequence = ''
for chr in sorted_chrs:
    sorted_positions = sorted(passed_snp_positions[chr])
    for pos in sorted_positions:
        if pos in amb_pos_counts[chr]:
            if amb_pos_counts[chr][pos] <= max_amb:
                sequence += ref_bases[chr][pos]
        else:
            sequence += ref_bases[chr][pos]

for i in range(0, len(sequence), 60):
    print(sequence[i:i+60])


# each sample
num_masked = 0
for sample in sample_list:
    print(">" + sample)
    sequence = ''
    for chr in sorted_chrs:
        sorted_positions = sorted(passed_snp_positions[chr])
        for pos in sorted_positions:
            if pos in amb_pos_counts[chr]:
                if amb_pos_counts[chr][pos] <= max_amb:
                    try:
                        if pos in alt_bases[sample][chr]:
                            sequence += alt_bases[sample][chr][pos]
                        else:
                            sequence += ref_bases[chr][pos]
                    except:
                        sequence += ref_bases[chr][pos]
                else:
                    num_masked += 1
            else:
                try:
                    if pos in alt_bases[sample][chr]:
                        sequence += alt_bases[sample][chr][pos]
                    else:
                        sequence += ref_bases[chr][pos]
                except:
                    sequence += ref_bases[chr][pos]
    for i in range(0, len(sequence), 60):
        print(sequence[i:i+60])

sys.stderr.write(str(num_masked)
                 + " sites filtered due to high missingness.\n")
