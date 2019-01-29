class fastq(analysis):
    def __init__(input, prefix=None, suffix=None, outdir='.', fasta=None,
                 gff=None, RAM=4, threads=1):
        analysis.__init__(input, prefix=None, suffix=None, outdir='.',
                          fasta=None, gff=None, RAM=4, threads=1)
        if len(self.input) == 1:
            self.paired = False
        elif len(self.input) == 2:
            self.paired = True
        else:
            raise ValueError('Input should only be one or two fastq files.')

    def bwa_align(self):
        """ perform bwa alignment
        :returns:
        """
        if not os.path.isfile(self.fasta+'.fai'):
            samtools_index_fa(self.fasta)
        if not (os.path.isfile(fa+'.bwt')):
            bwa_index_fa(fa)
        run(' '.join(['bwa mem', self.fasta, self.input, '| samtools view -S -b -u > ',
                      self.prefix+'.bam']))
        sorted_bam = sort_bam(self.prefix+'.bam')
        index_bam(sorted_bam)
        return sorted_bam

    def fastqc(self):
        """ quality control of raw or aligned BAMs using FastQc
        :return output directory
        """
        cmd = ' '.join(['fastqc -f fastq', self.input, '-o', self.outdir])
        run(cmd)
        print(' - FastQc finished.')
        return out_dir
