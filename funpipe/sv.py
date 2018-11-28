from . import utils.run

class sv(bam):
    def __init__(input, prefix=None, suffix=None, outdir='.', fasta=None,
                 gff=None, RAM=4, threads=1):
        bam.__init__(input, prefix=None, suffix=None, outdir='.', fasta=None,
                     gff=None, RAM=4, threads=1)

    def breakdancer(self):
        """ Detect structural variation using breakdancer
        :return
        """
        # create config files
        cfg_file = self.prefix+'.cfg'
        run('bam2cfg.pl -g -h'+self.input+'> '+prefix+'.cfg')
        # Detect chromosomal structural variants using breakdancer-max
        run('brakdancer-max -q 40 -r 20 -y 90 '+cfg_file)
        return cfg_file
