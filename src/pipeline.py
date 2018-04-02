from subprocess import check_call
import os
import contextlib

# from plumbum import local
# from plumbum.cmd import wget
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

def sort_bam(bam, out_dir='.', delete=False):
    ''' sort BAM using samtools
    :param bam: input bam
    :param out_dir: output directory, default '.'
    :param delete: delete original BAM
    '''
    bam_name = os.path.basename(bam)
    outfile = os.path.join(
        out_dir, os.path.splitext(bam_name)[0] + '.sorted.bam')
    run('samtools sort '+bam+' > '+outfile)
    if delete:
        rm(bam)
    return output

def samtools_index_fa(fa):
    ''' index a fasta file in place
    :param fa: fasta file
    '''
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
        samtools_index_fa(fa)
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
    fq1=prefix+'_1.fq.gz'
    fq2=prefix+'_2.fq.gz'
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
        '--vcf --changes --tracks --verbose > '+prefix+'.pilon.log 2>&1'])
    run(cmd)
    return cmd

def process_pilon_out(log, outdir, prefix):
    ''' process pilon output
        log: logfile
        outdir: output directory
    '''
    cmd = ' '.join(
         ['pilon_metrics', '-d', outdir, '-l', log, '--out_prefix', prefix])
    run(cmd)
    return cmd

def index_bam(bam):
    run('samtools index '+bam)
    return bam+'.bai'

def bwa_index_fa(fa):
    run('bwa index '+fa)
    return fa+'.fai'

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

def snpeff_db(
    gff3, dir, genome, config, prefix, ram, jar, ref_fa):
    ''' Create snpEff database
    gff3: gff file of gene annotation
    genome: name of the reference genome
    config: snpEff config files
    prefix: output Prefix
    ram: RAM in GB
    jar: snpEff jar
    ref_fa: reference fasta file
    '''
    snpeff_dir = os.path.dirname(jar)
    cmd = ' '.join(['sh snpeff_db.sh', dir, snpeff_dir, genome, ref_fa, gff3,
                    ram])
    run(cmd)
    return cmd

def bam_depth(bam, out_prefix, idx=False):
    ''' calculate bam depth
    :param bam: input to bam file
    :param out_prefix: output Prefix
    :param idx: output index
    '''
    outfile = out_prefix+'.depth.gz'
    cmd = 'samtools depth '+bam+' | bgzip > '+ outfile
    run(cmd)
    if idx:
        tabix(outfile, type='vcf')
    return outfile

def depth_per_window(pileup, out_prefix, faidx, window=5000):
    '''
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

def tabix(file, type=None):
    ''' Index tabix file
    :param file: input file
    :param type: file type, vcf
    '''
    cmd = 'tabix '+file
    if type:
        cmd += ' -p '+type
    run(cmd)

def fa2phylip(fa, out_prefix, jar):
    ''' transfer fasta file to phylip with java tool readSeq
    :param fa: fasta file
    :param jar: path to readseq.jar
    :param out_prefix:
    '''
    cmd = ' '.join(['java -cp', jar, 'run -f 12', fa])
    run(cmd)
    return

def ramxl(phylip, output, threads):
    ''' Run RAaML
    :param phylip: input phylip format file
    :param output: output file name
    :param threads: number of threads used for
    '''
    cmd = ' '.join([
        'raxmlHPC-PTHREADS-SSE3 -p 78960 -f a -x 12345 -N 1000 -m GTRCAT',
        '-T', str(threads), '-n', output, '-s', phylip])
    run(cmd)
    return output

def fasttree(fa, prefix):
    ''' Run FastTreeDP
    :param fa: fasta file
    :param prefix: output prefix
    '''
    cmd = ' '.join([
        'FastTreeDP -nt', fa, '>', out_prefix+'.nwk'
    ])
    return prefix+'.nwk'

def vcfsnpsToFa(vcf_list, out_prefix):
    out_file = out_prefix+'.fasta'
    cmd = ' '.join([
        'vcfSnpsToFasta.py', vcf_list, '>', out_prefix+'.fasta'
    ])
    run(cmd)
    return out_file

def filterGatkGenotypes(vcf, out_prefix):
    ''' filter Gatk output vcf
    :param vcf: input vcf file
    :param out_prefix: output prefix
    '''
    outfile = out_prefix+'_GQ50_AD08_DP10.vcf'
    cmd = ' '.join([
        'filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8',
        '--min_total_DP 10', vcf, '>', outfile
    ])
    run(cmd)
    return outfile


class picard:
    def __init__(self, jar='/seq/software/picard/1.853/bin/picard.jar', RAM=4):
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar
        ])

    def dict(self, fa, dict=None):
        ''' build fasta dictionary '''
        if dict is None:
            dict = 'tmp'
        cmd = ' '.join([self.cmd, "R="+fa, "O="+dict])
        run(cmd)

class gatk:
    def __init__(self, fa, jar='/xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar',
                 out_dir='.', RAM=4):
        ''' VCF sample QC
        :param vcf: vcf file
        :param fa: input Fasta
        :param jar: input jar file
        :param RAM: RAM usage
        :param out_dir: output directory
        '''
        self.out_dir = out_dir
        self.out = ''
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar, '-R', fa
        ])

    def variantEval(
        self, vcf, prefix, titv=True, samp=True, indel=True, multi=True):
        ''' VCF sample QC by different stratifications
        :param vcf: input vcf
        :param prefix: output prefix
        :param titv: use TiTv Evaluator
        :param indel: use InDel Evaluator
        :param multi: summarize multiallelic sites
        :param samp: stratify by samples
        '''
        out = os.path.join(self.out_dir, prefix+'.eval')
        cmd = ' '.join([self.cmd, '-T VariantEval', '--eval', vcf,
            '-o', out, '-noEV -noST -EV CountVariants'])
        if titv:
            cmd += ' -EV TiTvVariantEvaluator'
        if samp:
            cmd += ' -ST Sample'
        if indel:
            cmd += ' -EV IndelSummary'
        if multi:
            cmd += ' -EV MultiallelicSummary'
        run(cmd)
        return out

# def get_ref(ftp, md5, dir='.'):
#     """ download reference files from NCBI and perform md5sumcheck to files
#     :param ftp: ftp URL
#     :param md5: file that contains md5checksum results for
#     :param dir: destination directory
#     """
#     wget[ftp, dir]()
#     if md5:
#         with open(md5, 'r') as tsv:
#             for line in tsv:
#                 checksum, file = line.strip().split()
#                 check_md5(file, checksum)


# def check_md5(file, checksum):
#     """ perform md5 check sum
#     :param file: file to check
#     :param checksum: known checksum value
#     """
#     cmd = local['md5sum', file]
#     if :
#         return True
#     else:
#         return False
