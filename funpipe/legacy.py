import os
from math import ceil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

import sys
from funpipe.utils import run,rm
from funpipe.gt_pair import gt_pair
from funpipe.vcfheader import vcfheader
from funpipe.vcfrecord import vcfrecord
from funpipe.vcf import vcf


"""
Legacy methods for backward compatibility

"""
def pilon(fa, bam, prefix, ram, threads, jar):
    """ Run pilon commands

    Parameters
    ----------
        fa: :obj:`str` fasta file
        bam: :obj:`str` input bam path
        prefix: :obj:`str` output prefix
        ram: :obj:`int` input ram
        threads: :obj:`int` threads for pilon
        outdir: :obj:`str` output directory

    Returns
    -------


    """
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        '--genome', fa,
        '--frags', bam,
        '--output', prefix,
        '--threads', str(threads),
        '--vcf --changes --tracks --verbose > '+prefix+'.pilon.log 2>&1'])
    run(cmd)
    return cmd


def process_pilon_out(log, outdir, prefix):
    """ process pilon output
        log: logfile
        outdir: output directory
    """
    cmd = ' '.join(
         ['pilon_metrics', '-d', outdir, '-l', log, '--out_prefix', prefix])
    run(cmd)
    return cmd


def snpeff(invcf, outvcf, jar, config, genome, ram):
    """ run SNPEFF on a vcf
    invcf: input vcf
    outvcf: output vcf
    jar: snpeff jar
    genome: tag of genome name
    ram: memory in GB
    config: configuration file
    """
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        'eff', '-v',
        '-c', config,
        '-onlyCoding False',
        '-i vcf -o vcf', genome, invcf, '>', outvcf])
    run(cmd)
    return cmd


def snpeff_db(gff3, dir, genome, config, prefix, ram, jar, ref_fa):
    """ Create snpEff database
    gff3: gff file of gene annotation
    genome: name of the reference genome
    config: snpEff config files
    prefix: output Prefix
    ram: RAM in GB
    jar: snpEff jar
    ref_fa: reference fasta file
    """
    snpeff_dir = os.path.dirname(jar)
    cmd = ' '.join(['sh snpeff_db.sh', dir, snpeff_dir, genome, ref_fa, gff3,
                    ram])
    run(cmd)
    return cmd


def tabix(file, type=None):
    """ Index tabix file
    :param file: input file
    :param type: file type, vcf
    """
    cmd = 'tabix '+file
    if type:
        cmd += ' -p '+type
    run(cmd)
    return file+'.tbi'


def filterGatkGenotypes(vcf, out_prefix):
    """ filter Gatk output vcf
    :param vcf: input vcf file
    :param out_prefix: output prefix
    """
    outfile = out_prefix+'_GQ50_AD08_DP10.vcf'
    cmd = ' '.join([
        'filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8',
        '--min_total_DP 10', vcf, '>', outfile
    ])
    run(cmd)
    return outfile


def filter_variants(invcf, outvcf, min_GQ=50, AD=0.8, DP=10):
    """ apply variant filtering using GQ, AD and DP
    :param invcf: input vcf
    :param outvcf: output vcf
    :param min_GQ: minimum GQ cutoff
    :param AD: allelic depth cutoff
    :param DP: depth cutoff
    """
    cmd = ' '.join(['filterGatkGenotypes.py', '--min_GQ', str(min_GQ),
                    '--min_percent_alt_in_AD', str(AD),
                    '--min_total_DP', str(DP), invcf, '>', outvcf])
    run(cmd)
    return outvcf
