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

from scipy.stats import binom_test

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
        
        self.line3="20\t1234567\tmicrosat1\tGTCT\tG,GTACT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3"
        
        self.line4="chr1\t6\t.\tA\tG\t1659\tPASS\tDP=44;TD=77;BQ=38;MQ=38;QD=37;\
            BC=44,0,0,0;QP=100,0,0,0;PC=55;IC=0;DC=0;XC=1;AC=0;AF=0\tGT\t0/0"
        
        self.line5="chr1\t#\t#\t#\t#\t#\t#\tDP=44;TD=77;BQ=38;MQ=38;QD=37;\
            BC=44,0,0,0;QP=50,10,20,20;PC=55;IC=0;DC=0;XC=1;AC=0;AF=0\tGT\t#"
        
    def testPass(self):
        rec0 = vcfrecord(self.line0)
        rec1 = vcfrecord(self.line1)
        rec2 = vcfrecord(self.line2)
        self.assertTrue(  rec0.is_passing("GATK") )
        self.assertTrue( rec1.is_passing(".") )
        self.assertTrue( rec2.is_passing(".")==False )
    
    def testGetVarType(self):
        rec0 = vcfrecord(self.line0)
        rec1 = vcfrecord(self.line1)
        self.assertTrue(  rec0.get_variant_type(".", "0/1")=="SNP" )
        self.assertTrue( rec1.get_variant_type(".", "1/1")=="SNP" )
        
    def testGetVarLen(self):
        rec0 = vcfrecord(self.line0)
        rec1 = vcfrecord(self.line1)
        self.assertTrue(  rec0.get_variant_length("0")==False )
      
    
    #TODO: tests_npeff_annot
        
    def testGetAF(self):
        rec2 = vcfrecord(self.line2)
        rec4 = vcfrecord(self.line4)
        rec1 = vcfrecord(self.line1)
        self.assertTrue( rec4.get_AF() == 0 )
        self.assertTrue( rec2.get_AF() == 0.017 )
        self.assertTrue( rec1.get_AF() == 0.5 )
        
    def testGetQP(self):
        rec4 = vcfrecord(self.line4)
        self.assertDictEqual( rec4.get_QP(), {'A':100,'C':0,'G':0,'T':0} )
        rec5 = vcfrecord(self.line5)
        self.assertDictEqual( rec5.get_QP(), {'A':50,'C':10,'G':20,'T':20} )
        
    def testGetMAF(self):
        rec4 = vcfrecord(self.line4)
        self.assertTrue( rec4.get_MAF_from_QP() == 0 )
        rec5 = vcfrecord(self.line5)
        self.assertTrue( rec5.get_MAF_from_QP() == 0.5 )
    
    def testIsBia(self):
        rec0 = vcfrecord(self.line0)
        rec3 = vcfrecord(self.line3)
        self.assertTrue(  rec0.is_biallelic() == True  )
        self.assertTrue(  rec3.is_biallelic() == False )
        
    def testCountAmbig(self):
        rec0 = vcfrecord(self.line0)
        self.assertTrue(  rec0.count_ambig_genotypes() == 0 )
      
    def testGetGTProfile(self):
        rec0 = vcfrecord(self.line0)
        rec1 = vcfrecord(self.line1)
        rec2 = vcfrecord(self.line2)
        self.assertListEqual(  rec0.get_genotype_profile(),['0|0','0|0','0/1'] )
        self.assertListEqual(  rec1.get_genotype_profile(),['0|0','1|0','1/1'] )
        self.assertListEqual(  rec2.get_genotype_profile(),['0|0','0|1','0/0'] )
        
           

# test vcf
class TestVcf(unittest.TestCase):
    def setUp(self):
        self.test = 'test.vcf'
        self.test2 = 'test2.vcf'
        self.merged = 'merged.vcf'
        self.complex = 'phy_test.vcf'
        for file in ['output.dos.tsv','output.site_info.tsv','output.sample_info.tsv',
                    'merged.bed','merged.fam','merged.bim','merged.log',
                    'complex.bed','complex.fam','complex.bim','complex.log']:
            
            if os.path.exists(file):
                run('rm '+file)
                
                
    def test_n_samples(self):
        vcf_complex = vcf( self.complex )
        vcf_test = vcf(self.test)
        self.assertTrue( vcf_complex.n_samples == 4 )
        self.assertTrue(  vcf_test.n_samples == 1  )
    
    def testIndex(self):
        vcf_complex = vcf( self.complex )
        self.assertTrue( vcf_complex.get_sample_index( 'A' ) == 0 )
        self.assertTrue( vcf_complex.get_sample_index('B' ) == 1 )
        self.assertTrue( vcf_complex.get_sample_index('C' ) == 2 )
        self.assertTrue( vcf_complex.get_sample_index('D' ) == 3 )
        
    def testDos( self ):
        vcf_complex = vcf( self.complex )
        vcf_complex.cal_dos()
        self.assertTrue( os.path.exists( 'output.dos.tsv' ) )
        
    # site level test
    def testSiteInfo( self ):
        vcf_complex = vcf( self.complex )
        vcf_complex.get_site_info(info=['AF1'])
        self.assertTrue( os.path.exists( 'output.site_info.tsv' ) )
        self.assertTrue( vcf_complex.has_info( 'AF1', level='site' ) )
        
    # sample level test
    def testSampleInfo(self):
        vcf_complex = vcf( self.complex )
        vcf_complex.get_sample_info(info=['GT','SP'])
        self.assertTrue( os.path.exists( 'output.sample_info.tsv' ) )
        self.assertTrue( vcf_complex.has_info( 'GT', level='sample' ) )
        self.assertTrue( vcf_complex.has_info( 'SP', level='sample' ) )
        
    #test plink
    def testPlink_0(self):
        vcf_merged  = vcf( self.merged,prefix='merged'  )
        vcf_merged.get_plink()
        self.assertTrue( os.path.exists( 'merged.bed' ) )
        self.assertTrue( os.path.exists( 'merged.fam' ) )
        self.assertTrue( os.path.exists( 'merged.bim' ) )
        
    def testPlink_1(self):
        vcf_complex  = vcf( self.complex,prefix='complex'  )
        vcf_complex.get_plink()
        self.assertTrue( os.path.exists( 'complex.bed' ) )
        self.assertTrue( os.path.exists( 'complex.fam' ) )
        self.assertTrue( os.path.exists( 'complex.bim' ) )

#test GATK

# class TestGATK(unittest.TestCase):
#     def setUp(self):
#         self.test = 'test.vcf'
#         self.test2 = 'test2.vcf'
#         self.comp = 'test.comp.vcf'
      
       
#     def testCombineVar(self):
#         gatk_task = gatk('test.fa', prefix='combineVar')
#         vcf_dict = {'test':self.test,'test2':self.test2}
#         combined_vcf = gatk_task.combine_var(vcf_dict, option = 'UNSORTED', priority=None)
#         self.assertTrue(
#             filecmp.cmp(combined_vcf, 'merged.vcf')
#         )
        
       
#     def testVarEval(self):
#         gatk_task = gatk('test.fa', prefix='varEval')
#         eval_out = gatk_task.variant_eval( self.test )
#         self.assertTrue( os.path.exists( eval_out )  )
#         self.assertTrue( eval_out ==  './varEval.eval' )
        
      
#     def testSelectVar(self):
#         gatk_task = gatk('test.fa', prefix='selectVar')    
       


""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
