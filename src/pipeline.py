from subprocess import check_call
import os
import contextlib
# import configparser

@contextlib.contextmanager
def cd(dir):
    original_path = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(original_path)

def run(cmd):
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

def bam2fqs(bam, prefix, ram, picard):
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

def pilon(fa, bam, prefix, ram, threads, jar):
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

def update_snpeff_cfg(config, genome):
    """
    config: snpeff's configuration file
    genome: name of genome
    """

    return

def snpeff(invcf, outvcf, jar, config, genome, ram):
    ''' run SNPEFF on a vcf
    invcf: input vcf
    outvcf: output vcf
    jar: snpeff jar
    genome: tag of genome name
    ram: memory in GB
    config: configuration file
    '''
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        'eff', '-v',
        '-c', config,
        '-onlyCoding False',
        '-i vcf -o vcf', genome, invcf, '>', outvcf])
    run(cmd)
    return cmd

# def snpeff_db(gff3, genome, config, prefix, ram, jar):
#     ''' create snpEff database
#     gff: gff file of gene annotation
#     genome: name of the reference genome
#     config: snpEff config files
#     prefix: output Prefix
#     ram: RAM in GB
#     jar: snpEff jar
#     '''
#     os.mkdir(genome)
#     new_config = update_snpeff_cfg(config, genome)
#     cmd = ' '.join([
#         'java -Xmx4'+ram+'g',
#         '-jar', jar,
#         'build -gff3',
#         '-config', config])
#     return cmd


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
