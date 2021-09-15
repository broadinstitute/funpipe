import os
import sys
sys.path.append('.')
from picard import picard
from utils import run
from bam import bam
# from plumbum import local
# from plumbum.cmd import wget
# import configparser

class fastq:
    """ fastq """
    def __init__(self,*names,is_paired=False):
        '''Constructor of fastq object.
        
        Parameters
        ----------
        names: tuple
            The path(s) of fastq file(s).
        is_paired: bool
            True if there is a pair of fastq files, else False.
            
        Attributes
        ----------
        names: tuple of strings
            The path(s) of fastq file(s).
        is_paired: bool
            True if there is a pair of fastq files, else False.
        qc_outdir: string
            The output directory of fastq quality control report.
        bam: funpipe.bam
            The bam object of aligned sequences.
            
        '''
        self.is_paired = is_paired
        if len(names) > 2:
            raise Exception('Sorry, the maximum number of fastq files is 2')
            
        for fq in names:
            if not os.path.exists(fq):
                raise Exception('Sorry, input fastq file '+ fq + ' does not exist')
                
        self.names = names
        self.qc_outdir = None
        self.bam = None
        
    def fastqc(self,out_dir):
        ''' Quality control of raw or aligned BAMs using FastQc.
        
        Parameters
        ----------
        out_dir: string
            The output directory of quality control.
            
        Returns
        -------
        funpipe.fastq
            An updated fastq object with quality control files generated.
            
        '''
        if self.is_paired:
            fq1 = self.names[0]
            fq2 = self.names[1]
            cmd = ' '.join(['fastqc -f fastq', fq1, fq2, '-o',out_dir])
        else:
            fq = self.names[0]
            cmd = ' '.join(['fastqc -f fastq', fq, '-o', out_dir])
            
        run(cmd)
        print(' - FastQc finished.')
        self.qc_outdir = out_dir
        
        return self
        
        
    def bwa_align(self,fa,prefix):
        ''' Burrows-Wheeler alignment for sequences in fastq file(s).
        
        Parameters
        ----------
        fa: funpipe.fasta
            Reference fasta object.
        prefix: string
            Output file prefix.
            
        Returns
        -------
        funpipe.fastq
            An updated fastq object with the bam object of aligned sequences generated, path = prefix.bam.
            
        '''
        if fa.samtools_index == None or (not os.path.exists(fa.samtools_index)):
            fa.samtools_index_fa()
        if fa.bwa_index == None or (not os.path.exists(fa.bwa_index)):
            fa.bwa_index_fa()
        
        if self.is_paired:
            fq1 = self.names[0]
            fq2 = self.names[1]
            cmd = ' '.join(['bwa mem', fa.path, fq1, fq2, '| samtools view -S -b -u > ',
                      prefix+'.bam'])
        else:
            fq = self.names[0]
            cmd = ' '.join(['bwa mem', fa.path, fq, '| samtools view -S -b -u > ',
                      prefix+'.bam'])
            
        run( cmd )
        print(' - Burrows-Wheeler alignment finished.')
        self.bam = bam(prefix+'.bam')
        
        return self

