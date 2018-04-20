#!/usr/bin/env python3

from os.path import dirname, basename, splitext, join
import sys
import argparse
from crimson import picard
import json

pstats = {
    'alignment_summary_metrics':
        ['TOTAL_READS', 'PCT_PF_READS_ALIGNED', 'PCT_CHIMERAS'],
    'gc_bias.summary_metrics':
        ['AT_DROPOUT', 'GC_DROPOUT'],
    'insert_size_metrics':
        ['MEDIAN_INSERT_SIZE', 'STANDARD_DEVIATION'],
    'wgs_metrics': ['MEAN_COVERAGE']
}

stats_list = [
    'TOTAL_READS', 'PCT_PF_READS_ALIGNED', 'PCT_CHIMERAS', 'AT_DROPOUT',
    'GC_DROPOUT', 'MEDIAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'MEAN_COVERAGE']


def get_picard_stat(file):
    all_metr = picard.parse(file)['metrics']['contents']
    if isinstance(all_metr, dict):
        return all_metr
    else:
        return all_metr[-1]


def parse_file_path(path):
    ''' parse file path '''
    fdir = dirname(path)
    base = basename(path)
    prefix = splitext(base)[0]
    suffix = splitext(base)[-1]
    return fdir, base, prefix, suffix


def output_stats(qc_stats, bam_qc_file):
    with open(bam_qc_file, 'w+') as bam_qc:
        bam_qc.write('\t'.join(['Sample'] + stats_list)+'\n')
        for sample in qc_stats:
            stats = [sample]
            for stat in stats_list:
                stats.append(str(qc_stats[sample][stat]))
            bam_qc.write('\t'.join(stats)+'\n')


def process_bam_files(bam_list_file, bam_qc_file):
    '''
    :param bam_list_file: a list of bam files
    '''
    qc_stats = {}
    with open(bam_list_file, 'r') as bam_list:
        for line in bam_list:
            sample, path = line.strip().split('\t')
            fdir, fname, prefix, suffix = parse_file_path(path)
            qc_stats[sample] = {}
            for suffix in pstats:
                all_metr = get_picard_stat(join(fdir, prefix+'.'+suffix))
                for stat in pstats[suffix]:
                    qc_stats[sample][stat] = all_metr[stat]
    output_stats(qc_stats, bam_qc_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse Picard metrics from broad GP platform')
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='BAM list')
    required.add_argument(
        '-o', '--output', help="Output file")

    # optional arguments
    parser.add_argument(
        '-d', '--outdir', default='.', help='Output Directory')

    args = parser.parse_args()
    process_bam_files(
        args.input, join(args.outdir, args.output))
