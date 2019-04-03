import unittest
import sys
import os
from subprocess import check_call
import pandas as pd
import filecmp

from funpipe.utils import run
from funpipe.gatk import gatk
from scripts.coverage_analysis import get_contig_sets, \
     cal_chr_percent, cal_subg_percent, chr_coverage
from scripts.patch_ref_contigs import patch_fasta, patch_gff_contig


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


class TestPatchRefContigs(unittest.TestCase):
    # def test_chr_map_from_desc(self):
    #     desc_line = ''
    #     self.assertEqual()

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
"""


class TestCoverageAnalysis(unittest.TestCase):

    def setUp(self):
        self.all_contigs = {'1_A', '2_A', '1_D', '2_D', '1_B', '2_B'}
        self.subgenome = ['_A', '_B', '_D']

    def test_get_contig_sets(self):
        contig_map = get_contig_sets(self.all_contigs, self.subgenome)
        expect_out = {'1_A': 1, '2_A': 2, '1_B': 3, '2_B': 4, '1_D': 5,
                      '2_D': 6}
        self.assertDictEqual(contig_map, expect_out)

    def test_cal_chr_percent(self):
        cov = pd.DataFrame({
            'chr': ['chr1', 'chr1', 'chr2', 'chr2'],
            'sample0': [0, 1, 0, 1],
            'sample1': [1, 0, 1, 0]
        }, columns=['chr', 'sample0', 'sample1'])

        exp_cov_df = pd.DataFrame({
            'contigs': ['chr1', 'chr2'],
            'sample0': [0.5, 0.5],
            'sample1': [0.5, 0.5]
        }, columns=['contigs', 'sample0', 'sample1'])
        self.assertTrue(
            exp_cov_df.equals(cal_chr_percent(cov))
        )

    def test_cal_subg_percent(self):
        subg = ['_A', '_B']
        cov = pd.DataFrame({
            'chr': ['chr1_A', 'chr1_B', 'chr2_A', 'chr2_B'],
            'sample0': [0.4, 0.6, 0.4, 0.6],
            'sample1': [0.6, 0.4, 0.6, 0.4]
        }, columns=['chr', 'sample0', 'sample1'])
        exp_pct_df = pd.DataFrame({
            'subg': ['_A', '_B'],
            'sample0': [0.4, 0.6],
            'sample1': [0.6, 0.4]
        }, columns=['subg', 'sample0', 'sample1'])
        self.assertTrue(
            exp_pct_df.equals(cal_subg_percent(cov, subg))
        )

    def test_chr_coverage(self):
        cov_df = pd.DataFrame({
            'chr': ['chr1_A', 'chr1_A', 'chr2_A', 'chr2_A'],
            'sample0': [0.4, 0, 0.4, 0],
            'sample1': [0.6, 0.4, 0.6, 0.4]
        }, columns=['chr', 'sample0', 'sample1'])

        exp_cov_df = pd.DataFrame({
            'chr': ['chr1_A', 'chr2_A'],
            'sample0': [0.5, 0.5],
            'sample1': [1.0, 1.0]
        }, columns=['chr', 'sample0', 'sample1'])
        self.assertTrue(
            exp_cov_df.equals(chr_coverage(cov_df))
        )

    def test_pct_aneuploidy(self):
        cov_df = pd.DataFrame({
            'chr': ['chr1_A', 'chr1_A', 'chr2_A', 'chr2_A'],
            'sample0': [0.1, 1.2, 2.6, 3],
            'sample1': [3.9, 4.6, 2.1, 1.4]
        }, columns=['chr', 'sample0', 'sample1'])

        exp_cov_df = pd.DataFrame({
            'sample': ['sample0', 'sample0', 'sample0', 'sample0', 'sample0',
                       'sample1', 'sample1', 'sample1', 'sample1', 'sample1'],
            'ploidy': [  0,   1,   2,   3,   4,   0,   1,   2,   3,   4],
            'chr1_A': [0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            'chr2_A': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0]
        }, columns=['sample', 'ploidy', 'chr1_A', 'chr2_A'])

        out_cov_df = pct_aneuploidy(cov_df, max_ploidy=4)
        self.assertTrue(exp_cov_df.equals(out_cov_df))


'''
# test combineVariants:
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar \
    --variant test.vcf --variant test2.vcf -o merged.vcf \
    -genotypeMergeOptions UNSORTED -R test.fa -T CombineVariants \
    --assumeIdenticalSamples
'''


if __name__ == '__main__':
    os.chdir('data')
    unittest.main()
