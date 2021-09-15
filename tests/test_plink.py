import unittest
import sys
sys.path.append('../funpipe')
import os
from subprocess import check_call
import numpy as np
import pandas as pd
import filecmp

from utils import run
from plink import plink
from vcf import vcf

class TestPlink(unittest.TestCase):
    def setUp(self):
        
        self.bfile = ''
        
        
    def testConstructor(self):
        plink_test = plink(self.bfile)
        self.assertTrue( plink_test.bfile == self.bfile )
        self.assertTrue(  plink_test.related == None )
        self.assertTrue( plink_test.assoc == None )
        self.assertTrue( plink_test.qc == None )
        
#     def testRelatedness(self):
#         plink_test = plink(self.bfile)
#         plink_test.relatedness()
#         self.assertTrue( os.path.exists( plink_test.related )  )
#         self.assertTrue(  plink_test.related == self.bfile+'.related.tsv' )
        
#     def testGwas(self):
#         plink_test = plink(self.bfile)
#         plink_test.gwas()
#         self.assertTrue( os.path.exists( plink_test.assoc )  )
#         self.assertTrue(  plink_test.assoc == self.bfile + '.gemma.assoc.tsv' )
        
        
    def testImportPheno(self):
        plink_test = plink(self.bfile)
        pheno_list = []
        pheno_name = ''
        plink_test.import_pheno( pheno_list, pheno_name )
        temp = pd.read_csv( plink_test.fam, sep="\t")
        self.assertListEqual(  pheno_list, list(temp[pheno_name])   )
        
""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()       
        