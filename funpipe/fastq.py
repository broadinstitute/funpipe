from .picard import picard
from .utils import run
from .bam import bam
# from plumbum import local
# from plumbum.cmd import wget
# import configparser


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
        bam
            bam object containing indexed bam and sorted bam.
            
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
        out_dir = os.path.dirname(prefix)
        sorted_bam = BM.sort_bam( out_dir )
        indexed_bam = BM.index_bam()
        
        return BM

