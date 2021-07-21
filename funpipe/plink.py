"""
PLINK
=====
"""

class plink:
    def __init__(self, prefix):
        """constructor of plink object
        
        Parameters
        ----------
        arg1: string
            prefix: the prefix of bfile
        """
        self._bfile = prefix
        self._bed = prefix + '.bed'
        self._fam = prefix + '.fam'
        self._bim = prefix + '.bim'
        self._related = None
        self._assoc = None

    def relatedness(self):
        """ Calculate genetic relatedness matrix
        
        Returns
        -------
        string
            output relatedness file
            
        """
        self._related = self._bfile+'.related.tsv'
        run(" ".join(['gemma', '-bfile '+self._bfile,
                      '-gk -o '+self._related]))
        return self._related

    def gwas(self, lmm=4):
        """ Fit a LMM for a univariate model with GEMMA
        
        Parameters
        ----------
        arg1: int
            lmm: linear mixed model
            
        Returns
        -------
        string
            output gwas file
        """
        self._assoc = self._bfile + '.gemma.assoc.tsv'
        run(" ".join([
            'gemma -bfile '+self._bfile,
            '-k '+self._related,
            '-lmm '+str(lmm),
            '-o '+self._assoc]))
        return self._assoc

    def import_pheno():
        """ import phenotypes """

    def gwas_filter(self, ind=0.1, miss=0.1, maf=0.05):
        """ Filter genotype and sample level missingness and AF

        Parameters
        ----------
        arg1: float
            ind: individual level missingness, default = 0.1
        arg2: float
            miss: site level missingness, default = 0.1
        arg3: float
            maf: minor allele frequency cutoff, default = 0.05
            
        Returns
        -------
        string
            output filtered file

        """
        out_bfile = self._bfile+'.qc'
        run(" ".join([
            'plink --bfile', self._bfile,
            '--recode --out '+out_bfile,
            '--maf', maf,
            '--geno', miss,
            '--mind', ind
        ]))
        return out_bfile

    def variant_qc(self):
        #         run("plink --bfile --test-missing --allow-extra-chr")
        # plink --bfile ${PREFIX} --cluster missing --allow-extra-chr
        # plink --bfile ${PREFIX} --missing --allow-extra-chr
        return self

    def sample_qc(self):
        return self
