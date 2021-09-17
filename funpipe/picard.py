import os
import sys
from funpipe.utils import run
from funpipe.utils import rm

class picard:
    def __init__(self, jar='/opt/picard-tools/picard.jar', RAM=4):
        '''

        Parameters
        ----------
        jar: string
            The path to picard.jar, default = '/opt/picard-tools/picard.jar'.
        RAM: int
            Maximum RAM usage in gigabytes, default = 4.
            
        '''
        if not os.path.exists(jar):
            raise Exception('Sorry, jar file of picard tools does not exist')
            
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar
        ])


    def dict(self, fa, dictionary=None):
        ''' Build fasta dictionary.
        
        Parameters
        ----------
        fa: string
            The path to fasta file.
        dictionary: string
            The path to existing dictionary file, default = None.
            
        Returns
        -------
        string
            The path to existing dictionary file
            
        '''
        
        
        if dictionary is None:
            dictionary = os.path.splitext(fa)[0]+'.dict'
        
        if os.path.exists(dictionary):
            run('rm '+ dictionary)
            
            
        cmd = ' '.join([self.cmd, 'CreateSequenceDictionary', "R="+fa,
                       "O="+dictionary])
        run(cmd)
        
        return dictionary

    def bam2fqs(self, bam, prefix, is_paired = True ):
        ''' Realign BAM file to a reference genome, convert BAM to FASTQ.
        
        Parameters
        ----------
        bam: string
            The path to bam file.
        prefix: string
            The output prefix.
        is_paired: bool
            True if BAM is aligned from a pair of FASTQ files, else, False.
            
        Returns
        -------
        tuple
            The pair of FASTQ files or a single FASTQ file.
            
        '''
        out = None
        if is_paired:
            fq1 = prefix+'_1.fq.gz'
            fq2 = prefix+'_2.fq.gz'
            cmd = ' '.join([self.cmd, 'SamToFastq', 'I='+bam, 'F='+fq1, 'F2='+fq2])
            out = (fq1, fq2)
        else:
            fq = prefix+'.fq.gz'
            cmd = ' '.join([self.cmd, 'SamToFastq', 'I='+bam, 'F='+fq ])
            out = (fq)
            
        run(cmd)
        print("- bam2fqs Done.")
        
        return out
