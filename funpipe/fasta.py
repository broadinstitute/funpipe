from .picard import picard
from .utils import run


class fasta(analysis):
    def __init__(self, ):
        analysis.__init__()
        self.fa = fa
        return

    def bwa_index_fa(self):
        """Index a BWA file

        Examples
        --------

        >>> bwa_index_fa('~/examples/a.fa')  # doctest: +SKIP

        Parameters
        ----------
        fa : :obj:`str`
            fasta file

        Returns
        -------
        :obj:`str`
            Index file name
        """
        run('bwa index '+self.fa)
        return self.fa+'.fai'

    def samtools_index_fa(self):
        """Index a fasta file in place

        Examples
        --------

        >>> samtools_index_fa('~/examples/a.fa')  # doctest: +SKIP

        Parameters
        ----------
        fa : :obj:`str`
            fasta file

        Returns
        -------
            index file name
        """
        run('samtools faidx '+fa)
        return fa+'.fai'

    def index_fa(self):
        """Index fasta file with common genomic tools
        :param fa: fasta file
        :returns: None
        """
        samtools_index_fa(self.fa)
        bwa_index_fa(self.fa)
        pcd = picard()
        pcd.dict(self.fa)


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
