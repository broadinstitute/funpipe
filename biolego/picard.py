import os
from . import utils.run


class picard(analysis):
    def __init__(self, input, prefix=None, suffix=None, outdir='.', fasta=None,
                 gff=None, RAM=4, threads=1,
                 jar='/seq/software/picard/1.853/bin/picard.jar'):
        analysis.__init__(self, input, prefix=None, suffix=None, outdir='.',
                          fasta=None, gff=None, RAM=4, threads=1,)
        self.cmd = ' '.join(['java -Xmx'+str(RAM)+'g -jar', jar])

    def dict(self, dict=None):
        """ build fasta dictionary """
        if dict is None:
            dict = os.path.splitext(self.fasta)[0]+'.dict'
        cmd = ' '.join([self.cmd, 'CreateSequenceDictionary', "R="+self.fasta,
                       "O="+dict])
        run(cmd)

    def bam2fqs(self, paired=True):
        """ Realign BAM file to a reference genome
        :returns: tuple of fq file pairs
        """
        if paired:
            fq1 = self.prefix+'_1.fq.gz'
            fq2 = self.prefix+'_2.fq.gz'
            cmd = 'F='+fq1, 'F2='+fq2
            outs = (fq1, fq2)
        else:
            raise ValueError('Not yet support single end sequencing data.')
        cmd = ' '.join([self.cmd, 'SamToFastq', 'I='+self.input, output]) + cmd
        run(cmd)
        return outs
