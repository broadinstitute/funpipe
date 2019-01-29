import os
from .picard import picard as pcd
from .vcf import tabix
from .fasta import samtools_index_fa, bwa_index_fa
from .utils import run


class bam(analysis):
    def __init__(self, input, prefix, outdir='.', fasta, gff=None, RAM=4,
                 threads=1, suffix):
        """ initialize a BAM analysis

            :param fa: :obj:`str` fasta file
            :param bam: :obj:`str` input bam path
            :param prefix: :obj:`str` output prefix
            :param RAM: :obj:`int` input ram
            :param threads: :obj:`int` threads for pilon
            :param outdir: :obj:`str` output directory

        """
        analysis.__init__(input, prefix, outdir='.', fasta, gff=None, RAM=4,
                          threads=1)
        self.qc_stats = None

    def pilon(self, pilon_jar):
        """ Run pilon commands

        Parameters
        ----------
        pilon_jar

        Returns
        -------

        """
        cmd = ' '.join([
            'java -Xmx'+str(self.RAM)+'g',
            '-jar', self.pilon_jar,
            '--genome', self.fa,
            '--frags', self.input,
            '--output', self.prefix,
            '--threads', str(self.threads),
            '--vcf --changes --tracks --verbose > '+self.prefix
                +'.pilon.log 2>&1'])
        run(cmd)
        return cmd

    def index_bam(self):
        """
        Index a bam file

        Parameters
        ----------
        bam : :obj:`str`
            bam path
        Returns
        -------
        :obj:`str`
            bam index path
        """
        run('samtools index '+self.input)
        return self.input+'.bai'

    def sort_bam(self, tmp=None):
        """Sort BAM using samtools

        Parameters
        ----------
        bam: :obj:`str`
            input bam

        out_dir: :obj'`str`
            Output directory
        tmp: :obj'`str` temporary directory for

        Returns
        -------
        :obj:`str`
            Name of sorted bam
        """
        bam_name = os.path.basename(self.input)

        if tmp is None:
            tmp = self.prefix
        outfile = self.prefix + '.sorted.bam'
        run('samtools sort -T '+tmp+' '+bam+' > '+outfile)
        return outfile


    def bam_depth(self, idx=False):
        """ calculate bam depth
        :param idx: whether to index output mpileup file or not, default false
        """
        outfile = self.prefix+'.depth.gz'
        cmd = 'samtools depth '+self.input+' | bgzip > '+outfile
        run(cmd)
        if idx:
            tabix(outfile, type='vcf')
        return outfile

    def depth_per_window(self, pileup, faidx, window=5000):
        """ calculate depth per window
        :param faidx: fasta index
        :param pileup: pileup file from samtools
        :param window: window size in basepair
        """
        cmd = ' '.join([
            'dep_per_win.pl -m', pileup,
            '-p', self.prefix,
            '--window', str(window),
            '--faidx', faidx])
        run(cmd)
        return cmd

    def bam_sum(self):
        """ get read cateroties
        :return : output summary name
        """
        cmd = 'samtools flagstat ' + self.index + '> out_txt'
        run(cmd)
        return

    def clean_bam(self):
        """ clean up bams to keep only good quality alignment
        :param bam: input bam path
        :param out_prefix: output prefix
        """
        out_file = prefix+'.cleanup.bam'
        cmd = ' '.join(
            ['samtools view',
             '-bh -f 2 -F 1024 -F 4 -F 8 -F 512 -F 2048 -F 256 -q 30',
             self.input, '>', out_file]
        )
        run(cmd)
        return out_file
