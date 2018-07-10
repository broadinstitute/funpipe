import unittest
import sys
import os
from subprocess import check_call

from patch_ref_contigs import *
from post_process_coverage import get_contig_sets
from pipeline import gatk
import filecmp

# class TestMethod(unittest.TestCase):
#
#     def test_snpeff_db(self):
#
#         self.assertEqual()

os.chdir('data')


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


# class TestGetContigSets(unittest.TestCase):
#     def setUp(self):
#         self.all_contigs = {'1_A', '2_A', '1_D', '2_D', '1_B', '2_B'}
#         self.subgenome = ['_A', '_B', '_D']
#
#     def test_get_contig_sets(self):
#         contig_map = get_contig_sets(self.all_contigs, self.subgenome)
#         expect_out = {'1_A': 1, '2_A': 2, '1_B': 3, '2_B': 4, '1_D': 5,
#                       '2_D': 6}
#         self.assertDictEqual(contig_map, expect_out)
# test combineVariants: java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar --variant test.vcf --variant test2.vcf -o merged.vcf -genotypeMergeOptions UNSORTED -R test.fa -T CombineVariants --assumeIdenticalSamples




if __name__ == '__main__':
    unittest.main()
