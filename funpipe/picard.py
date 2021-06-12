import os
from .utils import run


class picard:
    def __init__(self, jar='/seq/software/picard/1.853/bin/picard.jar', RAM=4):
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar
        ])

    def dict(self, fa, dict=None):
        ''' build fasta dictionary '''
        if dict is None:
            dict = os.path.splitext(fa)[0]+'.dict'
        cmd = ' '.join([self.cmd, 'CreateSequenceDictionary', "R="+fa,
                       "O="+dict])
        run(cmd)

    def bam2fqs(self, bam, prefix):
        ''' Realign BAM file to a reference genome
        :param bam: bam file
        :param prefix: output prefix
        :parma fa: fasta file
        :returns: tuple of fq file pairs
        '''
        fq1 = prefix+'_1.fq.gz'
        fq2 = prefix+'_2.fq.gz'
        cmd = ' '.join([self.cmd, 'SamToFastq', 'I='+bam, 'F='+fq1, 'F2='+fq2])
        run(cmd)
        return(fq1, fq2)
