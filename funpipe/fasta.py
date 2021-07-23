from .picard import picard
from .utils import run
from .bam import bam
# from plumbum import local
# from plumbum.cmd import wget
# import configparser
"""
Fasta
=====
"""
class fasta:
    """fasta"""
    def __init__(self,fa_name):
        '''constructor of fasta object
        
        Parameters
        ----------
        fa_name: string
            the name of fasta file
        '''
        self.fa_name = fa_name
        self.picard_index = None
        self.bwa_index = None
        self.samtools_index =None
        
    def bwa_index_fa(self):
        ''' index a fasta file in place using bwa
        
        Returns
        -------
        string
            index file name
            
        '''
        run('bwa index '+self.fa_name)
        self.bwa_index = self.fa_name.split('.')[0]+'.bwt'
        return self.fa_name.split('.')[0]+'.bwt'


    def samtools_index_fa(self):
        ''' index a fasta file in place using samtools
        
        Returns
        -------
        string
            index file name
            
        '''
        run('samtools faidx '+self.fa_name)
        self.samtools_index = self.fa_name.split('.')[0]+'.fai'
        return self.fa_name.split('.')[0]+'.fai'


    def index_fa(self):
        ''' index fasta file with common genomic tools
        
        Returns
        -------
        string
            index file name
            
        '''
        samtools_index_fa()
        bwa_index_fa()
        pcd = picard()
        self.picard_index = pcd.dict( self.fa_name )


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

"""
FastQ
=====
"""
class fastq:
    """ fastq """
    def __init__(self,*names,is_paired=False):
        '''constructor of fastq object
        
        Parameters
        ----------
        names: tuple of strings
            name(s) of fastq file(s)
        is_paired: bool
            True if there is a pair of fastq files, else False.
            
        '''
        self.is_paired = is_paired
        self.names = names
        
        
    def fastqc(self,out_dir):
        ''' quality control of raw or aligned BAMs using FastQc
        
        Parameters
        ----------
        out_dir: string
            output directory
            
        Returns
        -------
        string
            output directory
            
        '''
        if self.is_paired:
            fq1 = self.names[0]
            fq2 = self.names[1]
            cmd = ' '.join(['fastqc -f fastq', fq1, fq2, '-o', out_dir])
        else:
            fq = self.names[0]
            cmd = ' '.join(['fastqc -f fastq', fq, '-o', out_dir])
            
        run(cmd)
        print(' - FastQc finished.')
        return out_dir
        
        
    def bwa_align(self,fa,prefix):
        ''' perform Burrows-Wheeler alignment
        
        Parameters
        ----------
        fa: string
            fasta file
        prefix: string
            output file prefix
            
        Returns
        -------
        string
            the name of sorted bam file
            
        '''
        if not os.path.isfile(fa+'.fai'):
            samtools_index_fa(fa)
        if not (os.path.isfile(fa+'.bwt')):
            bwa_index_fa(fa)
        
        if self.is_paired:
            fq1 = self.names[0]
            fq2 = self.names[1]
            cmd = ' '.join(['bwa mem', fa, fq1, fq2, '| samtools view -S -b -u > ',
                      prefix+'.bam'])
        else:
            fq = self.names[0]
            cmd = ' '.join(['bwa mem', fa, fq, '| samtools view -S -b -u > ',
                      prefix+'.bam'])
            
        run( cmd )
        BM = bam(prefix+'.bam')
        sorted_bam = BM.sort_bam(prefix+'.bam')
        BM.index_bam()
        return sorted_bam

