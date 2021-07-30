from .picard import picard
from .utils import run
from .bam import bam
# from plumbum import local
# from plumbum.cmd import wget
# import configparser


class fasta:
    """fasta"""
    def __init__(self,fa_name):
        '''constructor of fasta object
        
        Parameters
        ----------
        fa_name: string
            the name of fasta file
        '''
        self.fa_name = fa_name
        self.picard_index = None
        self.bwa_index = None
        self.samtools_index =None
        
    def bwa_index_fa(self):
        ''' index a fasta file in place using bwa
        
        Returns
        -------
        string
            index file name
            
        '''
        run('bwa index '+self.fa_name)
        self.bwa_index = self.fa_name.split('.')[0]+'.bwt'
        return self.fa_name.split('.')[0]+'.bwt'


    def samtools_index_fa(self):
        ''' index a fasta file in place using samtools
        
        Returns
        -------
        string
            index file name
            
        '''
        run('samtools faidx '+self.fa_name)
        self.samtools_index = self.fa_name.split('.')[0]+'.fai'
        return self.fa_name.split('.')[0]+'.fai'


    def index_fa(self):
        ''' index fasta file with common genomic tools
        
        Returns
        -------
        string
            index file name
            
        '''
        samtools_index_fa()
        bwa_index_fa()
        pcd = picard()
        self.picard_index = pcd.dict( self.fa_name )
        return self.picard_index


# def get_ref(ftp, md5, dir='.'):
#     """ download reference files from NCBI and perform md5sumcheck to files
#     :param ftp: ftp URL
#     :param md5: file that contains md5checksum results for
#     :param dir: destination directory
#     """
#     wget[ftp, dir]()
#     if md5:
#         with open(md5, 'r') as tsv:
#             for line in tsv:
#                 checksum, file = line.strip().split()
#                 check_md5(file, checksum)

