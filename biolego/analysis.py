import os

class analysis:
    """ bioanalysis class """
    def __init__(self, input, prefix='output', outdir='.', fasta=None, gff=None,
                 RAM=4, threads=1):
        """ initialize an analysis sample set

        Parameters
        ----------
        input:
            input file name
        prefix: prefix for this analysis
        outdir: output Directory
        RAM: RAM used for this analysis
        threads: number of threads used for this analysis

        """
        self.input = input
        self.prefix = os.path.join(outdir, prefix)
        self.outdir = outdir
        self.fasta = fasta
        self.gff = gff
        self.RAM = RAM
        self.threads = threads
