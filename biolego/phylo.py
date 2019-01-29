from . import utils.run


class phylo(analysis):
    def __init__(self, input, prefix='output', outdir='.', fasta=None, gff=None, RAM=4,
                 threads=1):
        analysis.__init__(self, input, prefix='output', outdir='.', fasta=None,
                          gff=None, RAM=4, threads=1)
        self.phylo_fa = None
        self.phylip = None
        self.ramxl = None
        self.free = None

    def vcf_to_fasta(self, max_amb=1000000):
        """ snp only vcf to fasta file
        :param max_amb: maximum number of samples with ambiguous calls for a
                        site to be included, recommended number of samples 10%
        """
        phylo_fa = self.prefix+'.fasta'
        cmd = ' '.join(['vcfSnpsToFasta.py --max_amb_samples', str(max_amb),
                        self.input, '>' phylo_fa])
        run(cmd)
        self.phylo_fa = phylo_fa
        return phylo_fa

    def fa2phylip(self, readseq_jar):
        """ transfer fasta file to phylip with java tool readSeq
        :param jar: path to readseq.jar
        """
        cmd = ' '.join(['java -cp', self.readseq_jar, 'run -f 12',
                        self.phylo_fa])
        self.phylip = self.prefix + '.phylip'
        run(cmd)
        return self.phylip

    def ramxl(self):
        """ Run RAMXL """
        self.ramxl = self.prefix + 'ramxl.txt'
        # phylip: input phylip format file
        cmd = ' '.join([
            'raxmlHPC-PTHREADS-SSE3 -p 78960 -f a -x 12345 -N 1000 -m GTRCAT',
            '-T', str(self.threads), '-n', self.ramxl, '-s', self.input])
        run(cmd)
        return self.ramxl

    def fasttree(self):
        """ Run FastTreeDP """
        self.free = self.prefix + '.nwk'
        cmd = ' '.join(['FastTreeDP -nt', self.phylo_fa, '>', self.free])
        run(cmd)
        return self.free
