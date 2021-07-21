import os
from .utils import run
"""
Picard
======
"""
class picard:
    def __init__(self, jar='/seq/software/picard/1.853/bin/picard.jar', RAM=4):
        '''
        
        Parameters
        ----------
        arg1: string
            jar: jar path of picard tools, default = '/seq/software/picard/1.853/bin/picard.jar'
        arg2: int
            RAM: RAM usage, default = 4
            
        '''
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar
        ])

    def dict(self, fa, dictionary=None):
        ''' build fasta dictionary
        
        Parameters
        ----------
        arg1: string
            fa: fasta file
        arg2: string
            dictionary: dictionary file, default = None
            
        Returns
        -------
        string
            dictionary file
            
        '''
        if dictionary is None:
            dictionary = os.path.splitext(fa)[0]+'.dict'
        cmd = ' '.join([self.cmd, 'CreateSequenceDictionary', "R="+fa,
                       "O="+dictionary])
        run(cmd)
        return dictionary

    def bam2fqs(self, bam, prefix):
        ''' Realign BAM file to a reference genome
        
        Parameters
        ----------
        arg1: string
            bam: bam file
        arg2: string
            prefix: output prefix
            
        Returns
        -------
        tuple
            tuple of fq file pairs
            
        '''
        fq1 = prefix+'_1.fq.gz'
        fq2 = prefix+'_2.fq.gz'
        cmd = ' '.join([self.cmd, 'SamToFastq', 'I='+bam, 'F='+fq1, 'F2='+fq2])
        run(cmd)
        return(fq1, fq2)
