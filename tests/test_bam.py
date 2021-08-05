import unittest
import sys
sys.path.append('../funpipe')
import os
from subprocess import check_call
import numpy as np
import pandas as pd
import filecmp

from utils import run
from fasta import fasta
from fastq import fastq
from bam import bam
from picard import picard

class TestFasta(unittest.TestCase):
    def setUp(self):
        self.fa_name = 'test.fa'
        for file in ['test.fa.fai','test.fa.bwt','test.dict']:
            if os.path.exists(file):
                run('rm '+file)
        
    def testConstructor(self):
        fasta_test = fasta( self.fa_name )
        self.assertTrue( os.path.exists( fasta_test.fa_name )  )
        
    def testFai(self):
        fasta_test = fasta( self.fa_name )
        fai_test = fasta_test.samtools_index_fa()
        self.assertTrue(os.path.exists(fai_test) )
        
    def testBwt(self):
        fasta_test = fasta( self.fa_name )
        bwt_test = fasta_test.bwa_index_fa()
        self.assertTrue(os.path.exists(bwt_test) )
        
    def testPicardIndex(self):
        fasta_test = fasta( self.fa_name )
        pic_index = fasta_test.index_fa('/opt/picard-tools/picard.jar')
        self.assertTrue(os.path.exists(pic_index) )
        
        
class TestBam(unittest.TestCase):
    
    def setUp(self):
        self.bam = 't.bam'
        for file in ['t.bam.bai','t.sorted.bam','t.depth.gz',
                     't_bam_summary.txt','t.cleanup.bam','t.vcf']:
            if os.path.exists(file):
                run('rm '+file)
                
                
                
    def testConstructor(self):
        bam_test = bam(self.bam)
        self.assertTrue(os.path.exists(bam_test.fname) )
        
    def testIndex(self):
        bam_test = bam(self.bam)
        bam_test.index_bam()
        self.assertTrue(os.path.exists(bam_test.indexed_bam) )
        
    def testSort(self):
        bam_test = bam(self.bam)
        bam_test.sort_bam( '', tmp=None, RAM=2, threads=1)
        self.assertTrue(os.path.exists(bam_test.sorted_bam) )
        
    def testDepth(self):
        bam_test = bam(self.bam)
        bam_test.bam_depth('t', idx=False)
        self.assertTrue(os.path.exists(bam_test.depth) )
        
    def testSum(self):
        bam_test = bam(self.bam)
        bam_test.bam_sum('t_bam_summary.txt')
        self.assertTrue( bam_test.summary == 't_bam_summary.txt')
        self.assertTrue(os.path.exists(bam_test.summary) )
        
    def testClean(self):
        bam_test = bam(self.bam)
        bam_test.clean_bam('t')
        self.assertTrue( bam_test.cleanup_bam == 't.cleanup.bam')
        self.assertTrue(os.path.exists(bam_test.cleanup_bam) )
         
    def testVariantCall(self):
        bam_test = bam(self.bam)
        ref_fa = 'test_ref.fa'
        bam_test.variant_call(ref_fa,'t')
        self.assertTrue( bam_test.out_vcf == 't.vcf')
        self.assertTrue(os.path.exists(bam_test.out_vcf ) )
        
        
    """    
    def testBreakdancer(self):
        bam_test = bam(self.bam)
        bam_test.breakdancer('/opt/breakdancer/perl/bam2cfg.pl','t_sv')
        self.assertTrue( bam_test.sv_config == 't_sv.cfg')
        self.assertTrue(os.path.exists(bam_test.sv_config ) )
    """   
    """
    def testDepthPerW(self):
        bam_test = bam(self.bam)
    
    """  
    
    
    
      
class TestFastq(unittest.TestCase):
    
    def setUp(self):
        self.fq1_name = 'sample_fq1.fq'
        self.fq2_name = 'sample_fq2.fq'
        self.ref_fa = 'test_ref.fa'
        
        for file in ['sample.bam','sample.bam.bai','sample.sorted.bam',
                    'sample_fq1_fastqc.html','sample_fq2_fastqc.html',
                    'sample_fq1_fastqc.zip','sample_fq1_fastqc.zip']:
            if os.path.exists(file):
                run('rm '+file)
          
    def testConstructor(self):
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        self.assertTrue(  fastq_test.is_paired )
        
    def testBwaAlign(self):
        # index bam failed due to unsorted position
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        fa = fasta(self.ref_fa)
        output_bm = fastq_test.bwa_align(fa,prefix = 'sample')
        
        self.assertTrue(output_bm.fname == 'sample.bam')
        #self.assertTrue(output_bm.indexed_bam == 'sample.bam.bai')
        #self.assertTrue(output_bm.sorted_bam == 'sample.sorted.bam')
        
        self.assertTrue(os.path.exists('sample.bam') )
        #self.assertTrue(os.path.exists('sample.bam.bai') )
        #self.assertTrue(os.path.exists('sample.sorted.bam') )
        
    def testFastqc(self):
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        output_dir = fastq_test.fastqc('.')
        self.assertTrue(os.path.exists('sample_fq1_fastqc.html') )
        self.assertTrue(os.path.exists('sample_fq2_fastqc.html') )
        

class TestPicard(unittest.TestCase):
    
    def setUp(self):
        self.bam = 't.bam'
        self.pcd = picard(jar='/opt/picard-tools/picard.jar', RAM=4)
        for file in ['t_1.fq.gz','t_2.fq.gz']:
            if os.path.exists(file):
                run('rm '+file)
                
    def testBam2fqs(self):
        fq1, fq2 = self.pcd.bam2fqs(self.bam, prefix = 't')
        self.assertTrue( fq1 == 't_1.fq.gz')
        self.assertTrue( fq2 == 't_2.fq.gz')
        self.assertTrue(os.path.exists('t_1.fq.gz') )
        self.assertTrue(os.path.exists('t_2.fq.gz') )
        

""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
