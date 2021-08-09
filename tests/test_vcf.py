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


# test vcf record


# test vcf




""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
