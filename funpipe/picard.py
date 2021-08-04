import os
import sys
sys.path.append('.')
from utils import run
from utils import rm
"""
Picard
======
"""
class picard:
    def __init__(self, jar='/opt/picard-tools/picard.jar', RAM=4):
        '''
        
        Parameters
        ----------
        jar: string
            jar path of picard tools, default = '/opt/picard-tools/picard.jar'
        RAM: int
            maximum RAM usage in gigabytes, default = 4
            
        '''
        if not os.path.exists(jar):
            raise Exception('Sorry, jar file of picard tools does not exist')
            
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar
        ])

    def dict(self, fa, dictionary=None):
        ''' build fasta dictionary
        
        Parameters
        ----------
        fa: string
            fasta file
        dictionary: string
            dictionary file, default = None
            
        Returns
        -------
        string
            dictionary file
            
        '''
        
        
        if dictionary is None:
            dictionary = os.path.splitext(fa)[0]+'.dict'
        
        if os.path.exists(dictionary):
            run('rm '+ dictionary)
            
            
        cmd = ' '.join([self.cmd, 'CreateSequenceDictionary', "R="+fa,
                       "O="+dictionary])
        run(cmd)
        
        return dictionary

    def bam2fqs(self, bam, prefix):
        ''' Realign BAM file to a reference genome
        
        Parameters
        ----------
        bam: string
            bam file
        prefix: string
            output prefix
            
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
