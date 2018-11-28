#!/usr/bin/env python3

from os.path import dirname, basename, splitext, join
import argparse
from typing import Any, Union

from crimson import picard
import pandas as pd
from glob import glob


stats = {
    'alignment_summary_metrics':
        ['TOTAL_READS', 'PCT_PF_READS_ALIGNED', 'PCT_CHIMERAS'],

    'wgs_metrics': ['MEAN_COVERAGE']
}

stats_list = [
    'TOTAL_READS', 'PCT_PF_READS_ALIGNED', 'PCT_CHIMERAS', 'MEAN_COVERAGE']


def realign_bam():
    fq1, fq2 = bam2fqs(args.bam, args.prefix, args.ram, args.picard_jar)
    realign_bam = bwa_align(args.fa, fq1, fq2, args.prefix)

def get_picard_stat(file):
    """

    :param file: input Picard file name
    :return:
    """
    all_metr = picard.parse(file)['metrics']['contents']
    if isinstance(all_metr, dict):
        return all_metr
    else:
        return all_metr[-1]


def get_ref_fa(analysis_file):
    """ get reference sequence from "analysis files"

    Examples
    --------

    >>> df = get_ref_fa('/seq/picard_aggregation/G143266/CL161/current/analysis_files.txt)
    :param analysis_file:
    """
    df = pd.read_csv(analysis_file, header=0, sep='\t')
    return df


def parse_gp_bam_path(path):
    """ parse on-prem BAM path from the genomic platform, to get directory,
    base name, prefix and suffix of the BAM. Prefix of the BAM is usually
    sample name.
    :param path: on-prem bam path from the genomic platform
    :return:
    """
    fdir = dirname(path)          # file directory
    base = basename(path)         # base name of the bam
    prefix = splitext(base)[0]    # prefix of the bam file, usually sample name
    suffix = splitext(base)[-1]   # suffix of the bam file
    return fdir, base, prefix, suffix


def output_stats(qc_stats, bam_qc_file):
    """ output qc statistics to a file
    :param qc_stats: an array that contains all BAM QC statistics
    :param bam_qc_file: file path to output all the metrics
    """
    with open(bam_qc_file, 'w+') as bam_qc:
        # output header to the file
        bam_qc.write('\t'.join(['Sample'] + stats_list)+'\n')
        # output QC stats for each sample from the array
        for sample in qc_stats:
            stats = [sample]
            for stat in stats_list:
                stats.append(str(qc_stats[sample][stat]))
            bam_qc.write('\t'.join(stats)+'\n')


def extract_picard_metrics(qc_path_tsv, bam_qc_file, is_gp_bam):
    """ process all bam files and get corresponding QC metrics and reference
    paths from Broad's Genomic Platform
    :param qc_path_tsv: a list of bam files
    :param bam_qc_file: file path to output all QC metrics
    :param is_gp_bam: input is a bam path from the Broad's GP
    """
    qc_stats = {}
    with open(qc_path_tsv) as qc_path:
        # process each bam record in the bam list
        for line in qc_path:
            sample, path = line.strip().split('\t')
            if is_gp_bam:
                fdir, fname, prefix, suffix = parse_gp_bam_path(path)
            qc_stats[sample] = {}
            for suffix in stats:
                stat_file = glob(join(path, sample+'.'+'*'+suffix))
                if len(stat_file) == 1:
                    all_metr: Union[dict, Any] = get_picard_stat(stat_file[0])
                    for stat in stats[suffix]:
                        qc_stats[sample][stat] = all_metr[stat]
                elif len(stat_file) > 1:
                        raise ValueError(sample+' has more than one '+suffix)
                else:
                    raise ValueError(suffix + ' not available.')
    output_stats(qc_stats, bam_qc_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            'Parse Picard metrics from a list of Picard output metrics. \n'
            'Note that this script assume file name is composed of '
            '<sample>.*.<suffix>'
        )
    )
    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', required=True, help='Input list of QC metrics')
    required.add_argument(
        '-o', '--output', help="Name of output file")

    # optional arguments
    parser.add_argument(
        '-d', '--outdir', default='.', help='Output Directory')
    parser.add_argument(
        '-b', '--is_gp_bam', action='store_true',
        help='whether input is a list of GP bam path'
    )
    args = parser.parse_args()
    extract_picard_metrics(
        args.input, join(args.outdir, args.output), args.is_gp_bam)
