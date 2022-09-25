import os
import sys
# unused
from funpipe.picard import picard as pcd
from funpipe.legacy import tabix
from funpipe.utils import run

class Bam:
    def __init__(self, path):
        '''
        
        Parameters
        ----------
        path: string
            The path to the BAM file.
            
        Attributes
        ----------
        path: string
            The path to BAM file.
        bam_index: string
            The path to index file.
        sorted_bam: funpipe.bam
            A bam object containing sorted BAM file.
        depth: string
            The path to depth file.
        pileup: string
            The path to pileup text file.
        depth_per_win: string
            The path to depth per window file.
        self.summary: string
            The path to summary text file.
        cleanup_bam: funpipe.bam
            A bam object containing cleanup BAM file.
        sv_config: string
            The path to config file containing structural variation detected by breakdancer.
        out_vcf: string
            The path to output VCF file by variant calling.

        Examples
        --------
        >>> from funpipe.bam import bam
        # sample = Bam('sample.bam')
        >>> bam = bam( 'sample.bam' )
        Sort BAM file:
        # sorted_bam = sample.sort_bam()
        >>> sorted_bam = bam.sort_bam().sorted_bam
        Clean BAM file:
        >>> cleanup_sorted_bam = sorted_bam.clean_bam().cleanup_bam
        Index BAM file:
        >>> cleanup_sorted_bam.index_bam()
        Compute the depth of alignment:
        >>> cleanup_sorted_bam.bam_depth( out_prefix = 'sample' )
        Detect structural variation:
        >>> cleanup_sorted_bam.breakdancer()
        Variant calling:
        >>> cleanup_sorted_bam.variant_calling( 'ref.fasta', 'sample')

        '''
        if os.path.exists(path):
            self.path = path
        else:
            raise Exception("Sorry, input bam file does not exist")
            
        self.bam_index = None
        self.sorted_bam = None
        self.depth = None
        self.pileup = None
        self.depth_per_win = None
        self.summary = None
        self.cleanup_bam = None
        self.sv_config = None
        self.out_vcf = None
    

    def sort_bam(self, out_dir=".", tmp=None, RAM=2, threads=1):
        ''' Sort BAM using samtools.
        
        Parameters
        ----------
        out_dir: string
            The output directory.
        tmp: string
            The temporary directory for the bam file, default is None.
        RAM: int
            Maximum RAM usage in gigabyte.
        threads: int
            Maximum number of threads used.
        
        Returns
        -------
        funpipe.bam
            An updated bam object with sorted bam object generated.
            
        '''
        bam_name = os.path.basename(self.path)
        prefix = os.path.join(out_dir, os.path.splitext(bam_name)[0])
        if tmp is None:
            tmp = prefix
        outfile = prefix + '.sorted.bam'
        run(' '.join(['samtools sort -T', tmp, '-m', str(RAM)+'G',
                      '-@', str(threads), '-o', outfile, self.path ]))
        # sorted_bam = outfile
        # return sorted_bam
        self.sorted_bam = bam(outfile)
        
        return self
    
    
    def index_bam(self):
        ''' Index BAM using samtools and generate ".bai" file.
        
        Notification
        ------------
        Index function of samtools might fail if BAM has unsorted positions.
        
        Returns
        -------
        funpipe.bam
            An updated bam object with indexed bam file generated.
            
        '''
        
        run('samtools index '+ self.path )
        # prefix name?
        self.bam_index = self.path +'.bai'
        
        return self


#    def bwa_align(fa, fq1, fq2, prefix):
#        ''' perform bwa alignment
#        :param fa: fasta file
#        :param fq1: fq file1
#        :param fq2: fq file2
#        :param prefix: output file prefix
#        :returns:
#        '''
#        if not os.path.isfile(fa+'.fai'):
#            samtools_index_fa(fa)
#        if not (os.path.isfile(fa+'.bwt')):
#            bwa_index_fa(fa)
#        run(' '.join(['bwa mem', fa, fq1, fq2, '| samtools view -S -b -u > ',
#                      prefix+'.bam']))
#        sorted_bam = sort_bam(prefix+'.bam')
#        index_bam(sorted_bam)
#        return sorted_bam

        

    def bam_depth(self, out_prefix, idx = False):
        '''Calculate bam depth.
        
        Parameters
        ----------
        out_prefix: string
            The output prefix, could include the output directory.
        idx: bool
            Whether to index output mpileup file or not, default is false.
            
        Returns
        -------
        funpipe.bam
            An updated bam object with depth file generated.
            
        '''
        outfile = out_prefix + '.depth.gz'
        cmd = 'samtools depth ' + self.path + ' | bgzip > ' + outfile
        run(cmd)
        if idx:
            tabix(outfile, type='vcf')
            
        self.depth = outfile 
        
        return self


#    def fastqc(bam, fq1, fq2, out_dir):
#        ''' quality control of raw or aligned BAMs using FastQc
#        :param bam: input bam path
#        :param out_dir: output directory
#        :param fq1: input fastq pair1
#        :param fq2: inpu  fastq pair2
#        :return output directory
#        '''
#        cmd = ' '.join(['fastqc -f fastq', fq1, fq2, '-o', out_dir])
#        run(cmd)
#        print(' - FastQc finished.')
#        return out_dir


    def create_pileup(self, fa, prefix, C = 50, reg = None, l = None, Q = 13, q = 0):
        """Produce pileup text format from an alignment
        
        Usage
        -----
        samtools mpileup [-EB] [-C capQcoef] [-r reg] [-f in.fa] [-l list] [-Q minBaseQ] [-q minMapQ] in.bam [in2.bam [...]]
        
        Parameters
        ----------
        fa: string
            Reference fasta file.
        out_prefix: string
            Prefix of output text file, '.txt' not included.
        C: int
            Coefficient for downgrading mapping quality for reads containing excessive mismatches.
            For BWA, default C = 50 is recommended.
        reg: string
            Only generate pileup in region. Requires the BAM files to be indexed.
            If not specified, default = None, pileup will be generated for all sites.
        l: string
            BED or position list file containing a list of regions or
            sites where pileup or BCF should be generated. Default = None.
        Q: int
            Minimum base quality for a base to be considered, default = 13.
        q: int
            Minimum mapping quality for an alignment to be used, default = 0.
        
        Returns
        -------
        funpipe.bam
            An updated bam object with pileup generated, out_prefix+'.txt'.
            
        """
        cmd = 'samtools mpileup -C ' + str(C)
        if reg != None:
            cmd = cmd + ' -r ' + reg
        cmd = cmd + ' -f ' + fa
        if l != None:
            cmd = cmd + ' -l ' + l
        cmd = cmd + ' -Q ' + str(Q) + ' -q ' + str(q) + ' -o ' + prefix + '.txt' + ' ' + self.path
        
        run(cmd)
        self.pileup = prefix + '.txt'
        
        return self
        
        
    def depth_per_window(self, out_prefix, faidx, window = 5000):
        '''Calculate depth per window.
        
        Parameters
        ----------
        pileup: string
            Pileup file from samtools.
        out_prefix: string
            The output prefix.
        faidx: string
            The path to fasta index file.
        window: int
            Window size in base pair.
            
        Returns
        -------
        funpipe.bam
            An updated bam object with depth per window file generated.
            
        '''
        if not os.path.exists( self.pileup ):
            raise Exception('Sorry, pileup file does not exist.')
        if not os.path.exists(faidx):
            raise Exception('Sorry, fasta index file does not exist.')
            
        cmd = ' '.join([
            'dep_per_win.pl -m', self.pileup,
            '-p', out_prefix,
            '--window', str(window),
            '--faidx', faidx])
        run(cmd)
        
        self.depth_per_win = out_prefix #assign depth per window
        
        return self


    def bam_sum(self,out_txt):
        '''Get read categories in the summary of BAM file.
        Count the number of alignments for each FLAG type.
        
        Parameters
        ----------
        out_txt: string
            The output file name,'.txt' must be included.
        
        Returns
        -------
        funpipe.bam
            An updated bam object with summary text generated.
            
        '''
        cmd = 'samtools flagstat ' + self.path + '>' + out_txt
        run(cmd)
        self.summary = out_txt
        
        return self


    def clean_bam(self, out_prefix):
        '''Clean up the BAM file to keep only good quality alignments.
        
        Parameters
        ----------
        out_prefix: string
            The output prefix.
            
        Returns
        -------
        funpipe.bam
            An updated bam object with cleaned bam object generated.
            
        '''
        out_file = out_prefix+'.cleanup.bam'
        cmd = ' '.join(
            ['samtools view',
             '-bh -f 2 -F 1024 -F 4 -F 8 -F 512 -F 2048 -F 256 -q 30',
             self.path, '>', out_file]
        )
        run(cmd)
        self.cleanup_bam = bam( out_file )
        
        return self
    
    
    def breakdancer(self, bd_path='/opt/breakdancer/', prefix='sv' ):
        ''' Detect structural variation using breakdancer.

        Parameters
        ----------
        bd_path: string
            The path to breakdancer, default = '/opt/breakdancer/'.
        prefix: string
            The output prefix.

        Returns
        -------
        funpipe.bam
            An updated bam object with structural variation detected by breakdancer.
            
        '''
        # create config files
        cfg_file = prefix+'.cfg'
        bam2cfg_path = os.path.join( bd_path , 'perl/bam2cfg.pl' )
        
        run( bam2cfg_path + ' -g -h '+ self.path +' > '+prefix+'.cfg')
        # Detect chromosomal structural variants using breakdancer-max
        run('breakdancer-max -q 40 -r 20 -y 90 '+cfg_file)
        self.sv_config = cfg_file
        
        return self
    
    def variant_call(self, fa, prefix):
        """Perform varaint calling from aligned bam file
        
        Parameters
        ----------
        fa: string
            Reference fasta file.    
        prefix: string
            The output prefix, without '.vcf' included.
            
        Returns
        -------
        funpipe.bam
            An updated bam object with variant calling executed, VCF file is generated.
        
        """
        #bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
        if not os.path.exists(fa):
            raise Exception('Sorry, reference fasta file for variant calling does not exist')
        
        output = prefix+'.vcf'
        
        cmd =' '.join(['bcftools mpileup','-f',fa,self.path,
                       '|','bcftools call -mv -o',prefix+'.vcf'])
        run(cmd)
        
        self.out_vcf = output
        
        return self
        