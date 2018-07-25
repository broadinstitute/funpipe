#!/usr/bin/env python3

from os.path import basename, splitext, join
import sys
import argparse
from pipeline import bam_depth, depth_per_window, sort_bam, cd, run, rm, \
    clean_bam, eprint


def run_ploidy(out_dir, bam, faidx, bam_sorted, clean_up):
    """ run ploidy analysis
    :param out_dir: output directory.
    :param bam: input BAM path.
    :param faidx: index of reference genome fasta.
    :param bam_sorted: whether the input BAM is sorted.
    :param chean_up: whether cleanup intemediate files.
    :return 1: workflow run successfully.
    """
    with cd(out_dir):
        base_prefix = splitext(basename(bam))[0]
        out_prefix = join(out_dir, base_prefix)
        bam_to_clean = bam
        if not bam_sorted:
            sorted_bam = sort_bam(bam, out_dir)
            bam_to_clean = sorted_bam
        cleaned_bam = clean_bam(bam_to_clean, out_prefix)
        pileup = bam_depth(cleaned_bam, out_prefix, idx=True)
        depth_per_window(pileup, base_prefix, faidx)
        if clean_up:  # clean up working directory
            if not bam_sorted:
                rm(sorted_bam)
            rm(cleaned_bam)
        eprint(" - Job finished.")
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate genomic ploidy for each bin from a BAM.'
        + 'BAMs are sorted and cleaned by keeping only high quality alignments'
        + '. Pileups were used to calculate dosage. Ploidy were calculated '
        + 'by normalizing average depth/bp against overall depth/bp.')
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-i', '--bam', required=True, help='Input BAM path'
    )
    required.add_argument(
        '--faidx', required=True, help='Fasta index'
    )

    # optional arguments
    parser.add_argument(
        '-o', '--out_dir', help="Output file", default='.'
    )
    parser.add_argument(
        '-n', '--bam_sorted', help='Do not sort input bam',
        action='store_true', default=False
    )
    parser.add_argument(
        '-c', '--cleanup', help='Cleanup intermediate files, mainly BAMs',
        action='store_true', default=False
    )

    args = parser.parse_args()
    run_ploidy(args.out_dir, args.bam, args.faidx, args.bam_sorted,
               args.cleanup)
