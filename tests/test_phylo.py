import unittest
import sys
sys.path.append('../funpipe')
import os
from subprocess import check_call
import numpy as np
import pandas as pd
import filecmp

from utils import run
from phylo import phylo
from vcf import vcf

class TestPhylo(unittest.TestCase):
    def setUp(self):
        self.invcf = 'phy_test.vcf'
        
    def testVcfSnp2Fasta(self):
        phylo_test = phylo(self.invcf)
        phylo_test.vcf_snp_to_fasta()
        
        self.assertTrue( os.path.exists( phylo_test.snp_fasta )  )
        self.assertTrue( phylo_test.snp_fasta == 'phy_test.snp.fasta' )
       
    def testFa2Phy(self):
        phylo_test = phylo(self.invcf)
        phylo_test.vcf_snp_to_fasta().fa2phylip()
        
        self.assertTrue( os.path.exists( phylo_test.phylip )  )
    
    def testRaxml(self):
        phylo_test = phylo(self.invcf)
        phylo_test.vcf_snp_to_fasta().fa2phylip().raxml()
        
        for out in phylo_test.raxml_result:
            self.assertTrue( os.path.exists( out ) )
            
    def testFastTree(self):
        phylo_test = phylo(self.invcf)
        phylo_test.vcf_snp_to_fasta().fa2phylip().fasttree()
        
        self.assertTrue( os.path.exists( phylo_test.fast_tree ) )
        
    def testPhyML(self):
        phylo_test = phylo(self.invcf)
        phylo_test.vcf_snp_to_fasta().fa2phylip().phyml()
        
        for out in phylo_test.phyml_result:
            self.assertTrue( os.path.exists( out ) )
            
            
            
"""main"""
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()       
        
        
        
        
        
        
