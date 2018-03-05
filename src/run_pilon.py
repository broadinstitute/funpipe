#!/usr/bin/env python3
import argparse
from subprocess import check_call
import os.path
import contextlib

@contextlib.contextmanager
def cd(dir):
    original_path = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(original_path)

def run(cmd):
    # print(cmd)
    check_call(cmd, shell=True)

def rm(file):
    run('rm '+file)

def sort_bam(bam, delete=True):
    '''
    bam: input bam
    delete: delete original bam
    '''
    output = bam[:-4] + '.sorted.bam'
    run('samtools sort '+bam+' > '+output)
    if delete:
        rm(bam)
    return output

def index_fa(fa):
    ''' index a fasta file in place '''
    run('samtools faidx '+fa)
    return fa+'.fai'

def bwa_align(fa, fq1, fq2, prefix):
    ''' perform bwa alignment
        fa: fasta file
        fq1: fq file1
        fq2: fq file2
        prefix: output file prefix
    '''
    if not os.path.isfile(fa+'.fai'):
        index_fa(fa)
    if not (os.path.isfile(fa+'.bwt')):
        bwa_index_fa(fa)
    run(' '.join(['bwa mem', fa, fq1, fq2, '| samtools view -S -b -u > ',
                  prefix+'.bam']))
    sorted_bam = sort_bam(prefix+'.bam')
    index_bam(sorted_bam)
    return sorted_bam

def bam2fqs(bam, prefix, ram=4,
            picard='/seq/software/picard-public/current/picard.jar'):
    ''' realign BAM file to a reference genome
        bam: bam file
        prefix: output prefix
        fa: fasta file
    '''
    fq1=prefix+'_1.fq'
    fq2=prefix+'_2.fq'
    cmd = ' '.join(['java -Xmx'+str(ram)+'g', '-jar', picard, 'SamToFastq',
                    'I='+bam, 'F='+fq1, 'F2='+fq2])
    run(cmd)
    return (fq1, fq2)

def pilon(fa, bam, prefix, ram, threads=1,
          jar = '/gsap/assembly_analysis/apps/prod/pilon/lib/pilon-1.12.jar'):
    '''
        fa: fasta file
        bam: bam
        prefix: output prefix
        ram: input ram
        threads: threads for pilon
        outdir: output directory
    '''
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        '--genome', fa,
        '--frags', bam,
        '--output', prefix,
        '--threads', str(threads),
        '--vcf --changes --tracks --verbose'])
    run(cmd)

def index_bam(bam):
    run('samtools index '+bam)

def bwa_index_fa(fa):
    run('bwa index '+fa)
    return fa+'.fai'

# def process_pilon_out(log, outdir):
#     ''' process pilon output
#         log: logfile
#         outdir: output directory
#     '''
#     script='/gsap/assembly_analysis/apps/prod/pilon/bin/pilon_metrics'
#     -l /path/pilon.logfile -d [/path/pilondirectory]
#     cmd = ' '.join([
#
#     ])

def main(bam, prefix, ram, fa, threads):
    fq1, fq2 = bam2fqs(bam, prefix, ram)
    realign_bam = bwa_align(fa, fq1, fq2, prefix)
    pilon(fa, realign_bam, prefix, ram, threads)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ram', default=16, type=int, help='RAM usage of input file')
    parser.add_argument(
        '--prefix', help='output prefix'
    )
    parser.add_argument(
        '--bam', help='input bam'
    )
    parser.add_argument(
        '--fa', help='reference genome fasta file to evaluate'
    )
    parser.add_argument(
        '--outdir', help='output directory'
    )
    parser.add_argument(
        '--threads', type=int, help='number of threads for pilon', default=4
    )
    args = parser.parse_args()
    with(cd(args.outdir)):
        main(args.bam, args.prefix, args.ram, args.fa, args.threads)
