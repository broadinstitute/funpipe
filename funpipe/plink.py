import os
import sys
sys.path.append('.')
from utils import run

class plink:
    def __init__(self, prefix):
        """constructor of plink object
        
        Parameters
        ----------
        prefix: string
            the prefix of bfile
        """
        self._bfile = prefix
        if os.path.exists(prefix + '.bed') and os.path.exists(prefix + '.fam') and os.path.exists(prefix + '.bim'):
            self._bed = prefix + '.bed'
            self._fam = prefix + '.fam'
            self._bim = prefix + '.bim'
        else:
            raise Exception('Sorry, missing at least one bfile for plink to work')
        
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
        lmm: int
            linear mixed model, 1 performs Wald test, 2 performs likelihood ratio test,
            3 performs score test, and 4 performs all the three tests, default = 4
            
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
        ind: float
            individual level missingness, default = 0.1
        miss: float
            site level missingness, default = 0.1
        maf: float
            minor allele frequency cutoff, default = 0.05
            
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
