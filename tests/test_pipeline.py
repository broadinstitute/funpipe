import unittest
import sys
import os
from subprocess import check_call

from patch_ref_contigs import *
import filecmp

# class TestMethod(unittest.TestCase):
#
#     def test_snpeff_db(self):
#
#         self.assertEqual()


class TestPatchRefContigs(unittest.TestCase):

    # def test_chr_map_from_desc(self):
    #     desc_line = ''
    #     self.assertEqual()
    os.chdir('data')

    def setUp(self):
        self.chr_map = {'NC_026745.1': 'chr1', 'NC_018792.1': 'chrMT'}

    def test_patch_fasta(self):
        patch_fasta('test.fa', 'outNoMt.fa', '_A', True)
        self.assertTrue(
            filecmp.cmp('outNoMt.fa', 'expNoMt.fa'),
            'Output and expect file not idential')

    def test_patch_gff_contig(self):
        patch_gff_contig(self.chr_map, 'test.gff', 'outNoMt.gff', '', True)
        self.assertTrue(
            filecmp.cmp('expNoMt.gff', 'outNoMt.gff'),
            'Output file not same as expected.')


if __name__ == '__main__':
    unittest.main()
