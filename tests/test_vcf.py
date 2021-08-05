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

# test gt_pair


#test vcf header


# test vcf record


# test vcf



""" main """
if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
