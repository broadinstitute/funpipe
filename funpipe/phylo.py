from funpipe.utils import run
"""
Phylogenetics
=============
"""
class phylo:
    def __init__(self, invcf ):
        """constructor of phylo object
        
        Parameters
        ----------
        invcf: string
            input vcf file
        """
        self.invcf = invcf
        
    def vcf_snp_to_fasta(self,prefix, max_amb=10):
        """ snp only vcf to fasta file

        Parameters
        ----------
        prefix: string
            output file prefix
        max_amb: int
            maximum number of samples with ambiguous calls for a site
            to be included, recommended number of samples 10%, use
            a very large number to disable this function 100000 (
            legacy options and will not be maintained.)
            
        Returns
        -------
        string
            output fasta file
        """
        cmd = ' '.join(['vcfSnpsToFasta.py --max_amb_samples', max_amb, self.invcf, '>',
                       prefix+'.fasta'])
        run(cmd)
        return prefix+'.fasta'


    def pairwise_snp_counts(self, fas, out_tsv):
        """ calculate pairwise snp overlap amongst a list of fasta using multiple
        alignment

        Parameters
        ----------
        fas: list
            a list of fasta files
        out_tsv: str
            output path of the pairwise SNP matrix
            
        Returns
        -------
        """
        return self

    def fa2phylip(self, fa, output, jar):
        ''' transfer fasta file to phylip with java tool readSeq
        
        Parameters
        ----------
        fa: string
            fasta file
        output: string
            output phylip file
        jar: string
            path to readseq.jar
        
        Returns
        -------
        string
            output phylip file
            
        '''
        cmd = ' '.join(['java -cp', jar, 'run -f 12', fa])
        run(cmd)
        return output


    def ramxl(self, phylip, output, threads):
        ''' Run RAxML:A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies
        
        Parameters
        ----------
        phylip: string
            input phylip format file
        output: string
            output file name
        threads: int
            number of threads used for
        
        Returns
        -------
        string
            output file name
            
        '''
        cmd = ' '.join([
            'raxmlHPC-PTHREADS-SSE3 -p 78960 -f a -x 12345 -N 1000 -m GTRCAT',
            '-T', str(threads), '-n', output, '-s', phylip])
        run(cmd)
        return output


    def fasttree(self, fa, prefix):
        ''' Run FastTreeDP: create phylogenetic trees from alignments of nucleotide or protein sequences
        
        Parameters
        ----------
        fa: string
            input fasta file
        prefix: string
            output prefix
            
        Returns
        -------
        string
            output phylogenetic tree file of newwick format
            
        '''
        cmd = ' '.join([
            'FastTreeDP -nt', fa, '>', prefix+'.nwk'
        ])
        run(cmd)
        return prefix+'.nwk'


#    def FastTreeDP(self, in_fa, out_prefix):
#        ''' perform FastaTreeDP analysis
#        
#        Parameters
#        ----------
#        in_fa: string
#            input fasta file
#        :param out_prefix: output file prefix
#        :returns nwk file
#        '''
#        out_nwk = out_prefix+'.nwk'
#        cmd = 'FastTreeDP -nt '+in_fa+' > ' + out_nwk
#        run(cmd)
#        return out_nwk

    # To do
    def phyml(self, phylip, output):
        """ Phylogenetic tree with parsimonies
        
        Parameters
        ----------
        phylip: string
            input phylip file
        output: string
            output file
        
        Returns
        -------
        string
            output file
            
        """
        cmd = ' '.join([
            'phyml -i', phylip
        ])
        return output
