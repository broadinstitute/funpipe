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
from gatk import gatk
from vcf import vcf
from gt_pair import gt_pair
from vcfheader import vcfheader
from vcfrecord import vcfrecord

# test gt_pair
class TestGT(unittest.TestCase):
    def setUp(self):
        self.gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        self.gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
    def testConstructor(self):
        gtP = gt_pair( self.gt1, self.gt2, False)
        self.assertListEqual( list(gtP.gt1)[0:5],list(self.gt1)[0:5] )
        self.assertListEqual( list(gtP.gt2)[0:5],list(self.gt2)[0:5] )
        
    def testGetTotal(self):
        gtP = gt_pair( self.gt1, self.gt2, False)
        gtP.get_n_total()
        self.assertTrue( gtP.n_total == 7 )
        
    def testGetShare(self):
        gtP = gt_pair( self.gt1, self.gt2, False)
        gtP.get_n_share()
        self.assertTrue( gtP.n_share == 2 )
        
    def testGetUnique(self):
        gtP = gt_pair( self.gt1, self.gt2, False)
        gtP.get_n_share()
        gtP.get_n_unique()
        self.assertTrue( gtP.n_unique == 5 )
        
        

#test vcf header
class TestVcfheader(unittest.TestCase):
    def setUp(self):
        self.vcf_header = vcfheader("test.vcf")
        
    def testGetSamples(self):
        self.assertListEqual( self.vcf_header.get_samples(), ["test.vcf"] )
        
    def testGetCaller(self):
        self.assertTrue( self.vcf_header.get_caller() == None )
        
    def testGetIndex(self):
        self.assertTrue( self.vcf_header.get_sample_index("test.vcf") == 0 )
        
    def testSnpeffStatus(self):
        self.assertTrue( self.vcf_header.get_snpeff_status() == False )
        
    def testGetContigs(self):
        self.assertListEqual(self.vcf_header.get_contigs(),["chr1","chr2"] )
    
    

# test vcf record
class TestVcfRecord(unittest.TestCase):
    def setUp(self):
        #3 examples from https://github.com/vcflib/vcflib/blob/master/samples/sample.vcf
        self.line0="19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3"
        
        self.line1="20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;\
        DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,."
        
        self.line2="20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.,."
        
    def testPass(self):
        rec0 = vcfrecord(self.line0)
        rec1 = vcfrecord(self.line1)
        rec2 =  vcfrecord(self.line2)
        self.assertTrue(  rec0.is_passing("GATK") )
        self.assertTrue( rec0.is_passing(".") )
        self.assertTrue( rec0.is_passing(".")==False )
        
    
    

# test vcf


#test GATK
"""
class TestGATK(unittest.TestCase):
    def setUp(self):
        self.eval = 'test.vcf'
        self.comp = 'test.comp.vcf'

    def testGenoConcord(self):
        gatk_task = gatk('test.fa', prefix='genoConcord')
        gatk_task.genotypeConcordance(self.comp, self.eval)
        self.assertTrue(
            filecmp.cmp('genoConcord.txt', 'genoConcord.exp.txt')
        )
        
    def testSelectVar(self):
        gatk_task = gatk('test.fa', prefix='selectVar')
        
    def testCombineVar(self):
        gatk_task = gatk('test.fa', prefix='combineVar')
    
    def testVarEval(self):
        gatk_task = gatk('test.fa', prefix='varEval')
        
"""        


""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
