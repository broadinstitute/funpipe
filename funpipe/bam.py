import os
import sys
sys.path.append('.')
from picard import picard as pcd
from vcf import tabix
from utils import run

class bam:
    """bam"""
    def __init__(self,filename):
        '''constructor of bam object
        
        Parameters
        ----------
        filename: string
            path to the bam file.
        '''
        if os.path.exists(filename):
            self.fname = filename
        else:
            raise Exception("Sorry, input bam file does not exist")
            
        self.indexed_bam = None
        self.sorted_bam = None
        self.depth = None
        self.depth_per_win = None
        self.summary = None
        self.cleanup_bam = None
        self.sv_config = None
        self.out_vcf = None
    
    def index_bam(self):
        ''' index BAM using samtools
        
        Returns
        -------
        string
            The name of indexed bam file
        '''
        run('samtools index '+ self.fname )
        
        self.indexed_bam = self.fname+'.bai'#assign indexed bam
        
        return self.indexed_bam


    def sort_bam(self,out_dir, tmp=None, RAM=2, threads=1):
        ''' sort BAM using samtools
        
        Parameters
        ----------
        out_dir: string
            output directory
        tmp: string
            temporary directory for the bam file, default is None
        RAM: int
            maximum RAM usage in gigabyte
        threads: int
            maximum number of threads used
        
        Returns
        -------
        string
            The name of sorted bam
            
        '''
        bam_name = os.path.basename(self.fname)
        prefix = os.path.join(out_dir, os.path.splitext(bam_name)[0])
        if tmp is None:
            tmp = prefix
        outfile = prefix + '.sorted.bam'
        run(' '.join(['samtools sort -T', tmp, '-m', str(RAM)+'G',
                      '-@', str(threads), '-o', outfile, self.fname ]))
        
        self.sorted_bam = outfile #assign sorted bam
        
        return self.sorted_bam


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

        

    def bam_depth(self, out_prefix, idx=False):
        ''' calculate bam depth
        
        Parameters
        ----------
        out_prefix: string
            output prefix, could include output path
        idx: bool
            whether to index output mpileup file or not, default is false
            
        Returns
        -------
        string
            The name of bam depth file
            
        '''
        outfile = out_prefix+'.depth.gz'
        cmd = 'samtools depth '+ self.fname +' | bgzip > '+outfile
        run(cmd)
        if idx:
            tabix(outfile, type='vcf')
            
        self.depth = outfile #assign bam depth
        
        return self.depth


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


    def depth_per_window(self, pileup, out_prefix, faidx, window=5000):
        ''' calculate depth per window
        
        Parameters
        ----------
        pileup: string
            pileup file from samtools
        out_prefix: string
            output prefix
        faidx: string
            fasta index file name
        window: int
            window size in basepair
            
        Returns
        -------
        string
            command line executed
        '''
        cmd = ' '.join([
            'dep_per_win.pl -m', pileup,
            '-p', out_prefix,
            '--window', str(window),
            '--faidx', faidx])
        run(cmd)
        
        self.depth_per_win = out_prefix #assign depth per window
        
        return self.depth_per_win


    def bam_sum(self,out_txt):
        ''' get read categories
        count the number of alignments for each FLAG type
        
        Parameters
        ----------
        out_txt: string
            output file name
        
        Returns
        -------
        string
            output summary name
        '''
        cmd = 'samtools flagstat ' + self.fname + '>' + out_txt
        run(cmd)
        self.summary = out_txt #assign bam summary
        return self.summary


    def clean_bam(self, out_prefix):
        ''' clean up bams to keep only good quality alignment
        
        Parameters
        ----------
        out_prefix: string
            output prefix
            
        Returns
        -------
        string
            output file name
            
        '''
        out_file = out_prefix+'.cleanup.bam'
        cmd = ' '.join(
            ['samtools view',
             '-bh -f 2 -F 1024 -F 4 -F 8 -F 512 -F 2048 -F 256 -q 30',
             self.fname, '>', out_file]
        )
        run(cmd)
        self.cleanup_bam = out_file #assign cleaned bam 
        
        return self.cleanup_bam
    
    
    def breakdancer(self, prefix):
        ''' Detect structural variation using breakdancer

        Parameters
        ----------
        prefix: string
            output prefix

        Returns
        -------
        string
            contig file
        '''
        # create config files
        cfg_file = prefix+'.cfg'
        run('bam2cfg.pl -g -h'+ self.fname +'> '+prefix+'.cfg')
        # Detect chromosomal structural variants using breakdancer-max
        run('brakdancer-max -q 40 -r 20 -y 90 '+cfg_file)
        self.sv_config = cfg_file
        
        return self.sv_config
    
    def variant_call(self,fa,prefix):
        """varaint calling from aligned bam file
        
        Parameters
        ----------
        fa: string
            reference fasta file
            
        prefix:
            output prefix, without '.vcf' included
            
        Returns
        -------
        string
            output vcf file
        
        Example
        -------
        bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
        
        """
        if not os.path.exists(fa):
            raise Exception('Sorry, reference fasta file for variant calling does not exist')
        
        output = prefix+'.vcf'
        
        cmd =' '.join(['bcftools mpileup','-f',fa,self.fname,
                       '|','bcftools call -mv -o',prefix+'.vcf'])
        run(cmd)
        
        self.out_vcf = output
        
        return self.out_vcf
        
        
        
        
        
        
        
        
    
    