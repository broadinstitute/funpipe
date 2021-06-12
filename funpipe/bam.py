import os
from .picard import picard as pcd
from .vcf import tabix
from .fasta import samtools_index_fa, bwa_index_fa
from .utils import run


def index_bam(bam):
    run('samtools index '+bam)
    return bam+'.bai'


def sort_bam(bam, out_dir, tmp=None, RAM=2, threads=1):
    ''' sort BAM using samtools
    :param bam: input bam
    :param out_dir: output directory
    :param tmp: temporary directory for
    :param RAM: maximum RAM usage in gigabyte
    :param t: maximum number of threads used
    :returns: name of sorted bam
    '''
    bam_name = os.path.basename(bam)
    prefix = os.path.join(out_dir, os.path.splitext(bam_name)[0])
    if tmp is None:
        tmp = prefix
    outfile = prefix + '.sorted.bam'
    run(' '.join(['samtools sort -T', tmp, '-m', str(RAM)+'G',
                  '-@', str(threads), '-o', outfile, bam]))
    return outfile


def bwa_align(fa, fq1, fq2, prefix):
    ''' perform bwa alignment
    :param fa: fasta file
    :param fq1: fq file1
    :param fq2: fq file2
    :param prefix: output file prefix
    :returns:
    '''
    if not os.path.isfile(fa+'.fai'):
        samtools_index_fa(fa)
    if not (os.path.isfile(fa+'.bwt')):
        bwa_index_fa(fa)
    run(' '.join(['bwa mem', fa, fq1, fq2, '| samtools view -S -b -u > ',
                  prefix+'.bam']))
    sorted_bam = sort_bam(prefix+'.bam')
    index_bam(sorted_bam)
    return sorted_bam


def bam_depth(bam, out_prefix, idx=False):
    ''' calculate bam depth
    :param bam: input to bam file
    :param out_prefix: output Prefix, could include output path
    :param idx: whether to index output mpileup file or not, default is false
    '''
    outfile = out_prefix+'.depth.gz'
    cmd = 'samtools depth '+bam+' | bgzip > '+outfile
    run(cmd)
    if idx:
        tabix(outfile, type='vcf')
    return outfile


def fastqc(bam, fq1, fq2, out_dir):
    ''' quality control of raw or aligned BAMs using FastQc
    :param bam: input bam path
    :param out_dir: output directory
    :param fq1: input fastq pair1
    :param fq2: inpu  fastq pair2
    :return output directory
    '''
    cmd = ' '.join(['fastqc -f fastq', fq1, fq2, '-o', out_dir])
    run(cmd)
    print(' - FastQc finished.')
    return out_dir


def depth_per_window(pileup, out_prefix, faidx, window=5000):
    ''' calculate depth per window
    :param pileup: pileup file from samtools
    :param window: window size in basepair
    '''
    cmd = ' '.join([
        'dep_per_win.pl -m', pileup,
        '-p', out_prefix,
        '--window', str(window),
        '--faidx', faidx])
    run(cmd)
    return cmd


def bam_sum(bam, out_txt):
    ''' get read cateroties
    :param bam: input bam path
    :param out_txt: output file name
    :return : output summary name
    '''
    cmd = 'samtools flagstat ' + bam + '> out_txt'
    run(cmd)
    return out_txt


def clean_bam(bam, out_prefix):
    ''' clean up bams to keep only good quality alignment
    :param bam: input bam path
    :param out_prefix: output prefix
    '''
    out_file = out_prefix+'.cleanup.bam'
    cmd = ' '.join(
        ['samtools view',
         '-bh -f 2 -F 1024 -F 4 -F 8 -F 512 -F 2048 -F 256 -q 30',
         bam, '>', out_file]
    )
    run(cmd)
    return out_file
