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
        self.path = 'test.fa'
        for file in ['test.fa.fai','test.fa.bwt','test.dict']:
            if os.path.exists(file):
                run('rm '+file)
        
    def testConstructor(self):
        fasta_test = fasta( self.path )
        self.assertTrue( os.path.exists( fasta_test.path )  )
        
    def testFai(self):
        fasta_test = fasta( self.path )
        fasta_test.samtools_index_fa()
        self.assertTrue(os.path.exists( fasta_test.samtools_index ) )
        
    def testBwt(self):
        fasta_test = fasta( self.path )
        fasta_test.bwa_index_fa()
        self.assertTrue(os.path.exists( fasta_test.bwa_index) )
        
    def testPicardIndex(self):
        fasta_test = fasta( self.path )
        fasta_test.index_fa('/opt/picard-tools/picard.jar')
        self.assertTrue(os.path.exists( fasta_test.picard_index) )
        
    def testFastaChain(self):
        fasta_test = fasta( self.path )
        fasta_test.samtools_index_fa().bwa_index_fa()
        self.assertTrue(os.path.exists(fasta_test.bwa_index) )
        self.assertTrue(os.path.exists(fasta_test.samtools_index) )
        
        
        
class TestBam(unittest.TestCase):
    
    def setUp(self):
        self.bam = 't.bam'
        self.bam_test = bam(self.bam)
        for file in ['t.sorted.bam.bai','t.sorted.bam','t.depth.gz',
                     't_bam_summary.txt','t.cleanup.bam','t.bam.bai',
                     't.vcf','t_pileup.txt','t.cleanup.sorted.bam','t.cleanup.sorted.bam.bai']:
            
            if os.path.exists(file):
                run('rm '+file)
                
                
                
    def testConstructor(self):
        self.assertTrue(os.path.exists(self.bam_test.path) )
        
    def testSort(self):
        self.bam_test.sort_bam( '', tmp=None, RAM=2, threads=1)
        self.assertTrue(os.path.exists( self.bam_test.sorted_bam.path ) )
        
    def testIndex(self):
        self.bam_test.index_bam()
        self.assertTrue( self.bam_test.bam_index == 't.bam.bai' )
        self.assertTrue(os.path.exists(self.bam_test.bam_index ) )
        
    def testDepth(self):
        self.bam_test.bam_depth('t', idx=False)
        self.assertTrue(os.path.exists(self.bam_test.depth) )
        
    def testSum(self):
        self.bam_test.bam_sum('t_bam_summary.txt')
        self.assertTrue( self.bam_test.summary == 't_bam_summary.txt')
        self.assertTrue( os.path.exists(self.bam_test.summary) )
        
    def testClean(self):
        self.bam_test.clean_bam('t')
        self.assertTrue( self.bam_test.cleanup_bam.path == 't.cleanup.bam')
        self.assertTrue(os.path.exists(self.bam_test.cleanup_bam.path) )
         
    def testVariantCall(self):
        ref_fa = 'test_ref.fa'
        self.bam_test.variant_call(ref_fa,'t')
        self.assertTrue( self.bam_test.out_vcf == 't.vcf')
        self.assertTrue(os.path.exists(self.bam_test.out_vcf ) )
        
    def testPileup(self):
        ref_fa = 'test_ref.fa'
        self.bam_test.create_pileup(fa=ref_fa,out_prefix='t_pileup',C=50,Q=13,q=0)
        self.assertTrue( self.bam_test.pileup == 't_pileup.txt')
        self.assertTrue(os.path.exists(self.bam_test.pileup) )
        
    """
    Error:
    ../../scripts/dep_per_win.pl --m t_pileup.txt --p t_depthPerW.txt --window 5000 --faidx test_ref.fa.fai
Experimental keys on scalar is now forbidden at ../../scripts/dep_per_win.pl line 31.
Type of arg 1 to keys must be hash or array (not hash element) at ../../scripts/dep_per_win.pl line 31, near "}) "
Execution of ../../scripts/dep_per_win.pl aborted due to compilation errors.
    

    def testDepthPerW(self):
        bam_test = bam(self.bam)
        
        ref_fa = 'test_ref.fa'
        fasta_test = fasta( ref_fa )
        fasta_test.samtools_index_fa()
        
        bam_test.create_pileup(fa=ref_fa,out_prefix='t_pileup',C=50,Q=13,q=0)
        depthPerW = bam_test.depth_per_window( 't_depthPerW',
                                              fasta_test.samtools_index, window=5000)
    """    
    
    """    
    def testBreakdancer(self):
        self.bam_test.breakdancer(prefix='t_sv')
        self.assertTrue( self.bam_test.sv_config == 't_sv.cfg')
        self.assertTrue(os.path.exists(self.bam_test.sv_config ) )
    """
    def testBAMChain_0(self):
        if os.path.exists('t.bam.bai'):
            run('rm ' + 't.bam.bai' )
        
        if os.path.exists('t.sorted.bam'):
            run('rm ' + 't.sorted.bam' )
            
        self.bam_test.index_bam().sort_bam( '', tmp=None, RAM=2, threads=1)
        self.assertTrue(os.path.exists(self.bam_test.bam_index ) )
        self.assertTrue(os.path.exists(self.bam_test.sorted_bam.path) )
        
    def testBAMChain_1(self):
        for file in ['t.depth.gz','t_bam_summary.txt','t.cleanup.bam']:
            if os.path.exists(file):
                run('rm ' + file)
                
        self.bam_test.bam_depth('t', idx=False).bam_sum('t_bam_summary.txt').clean_bam('t')
        self.assertTrue(os.path.exists(self.bam_test.depth) )
        self.assertTrue(os.path.exists(self.bam_test.summary) )
        self.assertTrue(os.path.exists(self.bam_test.cleanup_bam.path) )
        
    def testBAMChain_2(self):
        for file in ['t.sorted.bam','t.sorted.bam.bai']:
            if os.path.exists(file):
                run('rm ' + file)
                
        sorted_bam = self.bam_test.sort_bam( '', tmp=None, RAM=2, threads=1).sorted_bam
        sorted_bam.index_bam()
        self.assertTrue(os.path.exists(sorted_bam.path) )
        self.assertTrue(os.path.exists(sorted_bam.bam_index) )
        
    def testBAMChain_3(self):
        for file in ['t.cleanup.bam','t.cleanup.sorted.bam','t.cleanup.sorted.bam.bai']:
            if os.path.exists(file):
                run('rm ' + file )
                
        cleanup_bam = self.bam_test.clean_bam('t').cleanup_bam
        sorted_cleanup_bam = cleanup_bam.sort_bam(out_dir="",).sorted_bam
        sorted_cleanup_bam.index_bam()
        
        self.assertTrue(os.path.exists(cleanup_bam.path) )
        self.assertTrue(os.path.exists( sorted_cleanup_bam.path ) )
        self.assertTrue(os.path.exists(sorted_cleanup_bam.bam_index) )
        
        self.assertTrue(cleanup_bam.path == 't.cleanup.bam' )
        self.assertTrue( sorted_cleanup_bam.path == 't.cleanup.sorted.bam' )
        self.assertTrue( sorted_cleanup_bam.bam_index == 't.cleanup.sorted.bam.bai' )
         
    
      
class TestFastq(unittest.TestCase):
    
    def setUp(self):
        self.fq1_name = 'sample_fq1.fq'
        self.fq2_name = 'sample_fq2.fq'
        self.ref_fa = 'test_ref.fa'
        
        for file in ['sample.bam','sample.sorted.bam','sample.sorted.bam.bai',
                    'sample_fq1_fastqc.html','sample_fq2_fastqc.html',
                    'sample_fq1_fastqc.zip','sample_fq1_fastqc.zip']:
            if os.path.exists(file):
                run('rm '+file)
          
    def testConstructor(self):
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        self.assertTrue(  fastq_test.is_paired )
        
    def testBwaAlign(self):
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        fa = fasta(self.ref_fa)
        
        output_bm = fastq_test.bwa_align(fa,prefix = 'sample').bam
        sorted_bam = output_bm.sort_bam(out_dir="").sorted_bam
        sorted_bam.index_bam()
        
        self.assertTrue(output_bm.path == 'sample.bam')
        self.assertTrue(sorted_bam.path == 'sample.sorted.bam')
        self.assertTrue(sorted_bam.bam_index == 'sample.sorted.bam.bai')
        
        self.assertTrue(os.path.exists('sample.bam') )
        self.assertTrue(os.path.exists('sample.sorted.bam') )
        self.assertTrue(os.path.exists('sample.sorted.bam.bai') )
        
        
    def testFastqc(self):
        fastq_test = fastq(self.fq1_name,self.fq2_name,is_paired=True)
        fastq_test.fastqc('.')
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
