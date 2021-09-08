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


class TestBam(unittest.TestCase):
    def setUp(self):
        self.bam = 't.bam'
        self.bam_test = bam(self.bam)
        for file in ['t_pileup.txt']:
            
            if os.path.exists(file):
                run('rm '+file)
                
    def testDepthPerW(self):
        bam_test = bam(self.bam)
        
        ref_fa = 'test_ref.fa'
        fasta_test = fasta( ref_fa )
        fasta_test.samtools_index_fa()
        
        bam_test.create_pileup(fa=ref_fa,out_prefix='t_pileup',C=50,Q=13,q=0)
        depthPerW = bam_test.depth_per_window( 't_depthPerW',fasta_test.samtools_index, window=5000)
        
        
        
        
""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()