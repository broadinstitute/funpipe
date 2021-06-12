from .picard import picard
from .utils import run
# from plumbum import local
# from plumbum.cmd import wget
# import configparser


def bwa_index_fa(fa):
    run('bwa index '+fa)
    return fa+'.fai'


def samtools_index_fa(fa):
    ''' index a fasta file in place
    :param fa: fasta file
    :returns: index file name
    '''
    run('samtools faidx '+fa)
    return fa+'.fai'


def index_fa(fa):
    ''' index fasta file with common genomic tools
    :param fa: fasta file
    :returns: None
    '''
    samtools_index_fa(fa)
    bwa_index_fa(fa)
    pcd = picard()
    pcd.dict(fa)


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
