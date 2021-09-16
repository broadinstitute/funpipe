import os
import sys
import pandas as pd
from funpipe.utils import run

class plink:
    def __init__(self, prefix):
        """
        
        Parameters
        ----------
        prefix: string
            The prefix of bfile.
            
        Attributes
        ----------
        bed: string
            The path to bed file.
        fam: string
            The path to fam file.
        bim: string
            The path to bim file.
        related: string
            The path to genetic relatedness matrix.
        assoc: string
            The path to gwas result file.
        qc: string
            The path to quality control file.
        
        Examples
        --------
        >>> from funpipe.plink import plink
        >>> plink = plink( prefix = 'sample' )
        Perform GWAS:
        >>> plink.relatedness().gwas()

        """
        self._bfile = prefix
        if not os.path.exists(prefix + '.bed'):
            raise Exception('Sorry, missing bed file for plink to work')
        if not os.path.exists(prefix + '.fam'):
            raise Exception('Sorry, missing fam file for plink to work')
        if not os.path.exists(prefix + '.bim'):
            raise Exception('Sorry, missing bim file for plink to work')
        
        self._bed = prefix + '.bed'
        self._fam = prefix + '.fam'
        self._bim = prefix + '.bim'
        
        self._related = None
        self._assoc = None
        self._qc = None
        
    @property
    def related(self):
        return self._related
    @property
    def assoc(self):
        return self._assoc
    @property
    def qc(self):
        return self._qc
    @property
    def bfile(self):
        return self._bfile
    @property
    def bed(self):
        return self._bed
    @property
    def fam(self):
        return self._fam
    @property
    def bim(self):
        return self._bim
    
    def relatedness(self):
        """ Calculate genetic relatedness matrix.
        
        Returns
        -------
        funpipe.plink
            An updated plink object with relatedness file generated.
            
        """
        self._related = self._bfile+'.related.tsv'
        run(" ".join(['gemma', '-bfile '+self._bfile,
                      '-gk -o '+self._related]))
        
        return self

    def gwas(self, lmm=4):
        """ Fit a LMM for a univariate model with GEMMA.
        
        Parameters
        ----------
        lmm: int
            Linear mixed model, 1 performs Wald test, 2 performs likelihood ratio test,
            3 performs score test, and 4 performs all the three tests, default = 4.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with gwas result file generated.
            
        """
        self._assoc = self._bfile + '.gemma.assoc.tsv'
        run(" ".join([
            'gemma -bfile '+self._bfile,
            '-k '+self._related,
            '-lmm '+str(lmm),
            '-o '+self._assoc]))
        
        return self

    def import_pheno(self,phenotypes ):
        """ Import phenotypes
        
        Parameters
        ----------
        phenotypes: list
            The list of phenotypes appended to fam file.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with phenotype imported in fam file.
            
        """
        fam_pd = pd.read_csv(self._fam , sep = "\t", header = None)
        if len(fam_pd[ list(fam_pd.columns)[0] ]) != len(phenotypes):
            raise Exception("Unmatched length between imported phenotypes and fam file.")
        else:
            fam_pd[len(fam_pd.columns)] = phenotypes
        fam_pd.to_csv( self._bfile + '.fam', sep = "\t", index = None)
        self._fam = self._bfile + '.fam'
        
        return self

    def gwas_filter(self, ind=0.1, miss=0.1, maf=0.05):
        """ Filter genotype and sample level missingness and MAF.

        Parameters
        ----------
        ind: float
            Individual level missingness, default = 0.1.
        miss: float
            Site level missingness, default = 0.1.
        maf: float
            Minor allele frequency cutoff, default = 0.05.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with quality control file generated.

        """
        self._qc = self._bfile+'.qc'
        run(" ".join([
            'plink2 --bfile', self._bfile,
            '--recode --out '+ self._qc ,
            '--maf', maf,
            '--geno', miss,
            '--mind', ind
        ]))
        return self

#    def variant_qc(self):
        #         run("plink --bfile --test-missing --allow-extra-chr")
        # plink --bfile ${PREFIX} --cluster missing --allow-extra-chr
        # plink --bfile ${PREFIX} --missing --allow-extra-chr
#        return self

#    def sample_qc(self):
#        return self
import os
import sys
import pandas as pd
from funpipe.utils import run

class plink:
    def __init__(self, prefix):
        """
        
        Parameters
        ----------
        prefix: string
            The prefix of bfile.
            
        Attributes
        ----------
        bed: string
            The path to bed file.
        fam: string
            The path to fam file.
        bim: string
            The path to bim file.
        related: string
            The path to genetic relatedness matrix.
        assoc: string
            The path to gwas result file.
        qc: string
            The path to quality control file.
        
        Examples
        --------
        >>> from funpipe.plink import plink
        >>> plink = plink( prefix = 'sample' )
        Perform GWAS:
        >>> plink.relatedness().gwas()

        """
        self._bfile = prefix
        if not os.path.exists(prefix + '.bed'):
            raise Exception('Sorry, missing bed file for plink to work')
        if not os.path.exists(prefix + '.fam'):
            raise Exception('Sorry, missing fam file for plink to work')
        if not os.path.exists(prefix + '.bim'):
            raise Exception('Sorry, missing bim file for plink to work')
        
        self._bed = prefix + '.bed'
        self._fam = prefix + '.fam'
        self._bim = prefix + '.bim'
        
        self._related = None
        self._assoc = None
        self._qc = None
        
    @property
    def related(self):
        return self._related
    @property
    def assoc(self):
        return self._assoc
    @property
    def qc(self):
        return self._qc
    @property
    def bfile(self):
        return self._bfile
    @property
    def bed(self):
        return self._bed
    @property
    def fam(self):
        return self._fam
    @property
    def bim(self):
        return self._bim
    
    def relatedness(self):
        """ Calculate genetic relatedness matrix.
        
        Returns
        -------
        funpipe.plink
            An updated plink object with relatedness file generated.
            
        """
        self._related = self._bfile+'.related.tsv'
        run(" ".join(['gemma', '-bfile '+self._bfile,
                      '-gk -o '+self._related]))
        
        return self

    def gwas(self, lmm=4):
        """ Fit a LMM for a univariate model with GEMMA.
        
        Parameters
        ----------
        lmm: int
            Linear mixed model, 1 performs Wald test, 2 performs likelihood ratio test,
            3 performs score test, and 4 performs all the three tests, default = 4.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with gwas result file generated.
            
        """
        self._assoc = self._bfile + '.gemma.assoc.tsv'
        run(" ".join([
            'gemma -bfile '+self._bfile,
            '-k '+self._related,
            '-lmm '+str(lmm),
            '-o '+self._assoc]))
        
        return self

    def import_pheno(self,phenotypes ):
        """ Import phenotypes
        
        Parameters
        ----------
        phenotypes: list
            The list of phenotypes appended to fam file.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with phenotype imported in fam file.
            
        """
        fam_pd = pd.read_csv(self._fam , sep = "\t", header = None)
        if len(fam_pd[ list(fam_pd.columns)[0] ]) != len(phenotypes):
            raise Exception("Unmatched length between imported phenotypes and fam file.")
        else:
            fam_pd[len(fam_pd.columns)] = phenotypes
        fam_pd.to_csv( self._bfile + '.fam', sep = "\t", index = None)
        self._fam = self._bfile + '.fam'
        
        return self

    def gwas_filter(self, ind=0.1, miss=0.1, maf=0.05):
        """ Filter genotype and sample level missingness and MAF.

        Parameters
        ----------
        ind: float
            Individual level missingness, default = 0.1.
        miss: float
            Site level missingness, default = 0.1.
        maf: float
            Minor allele frequency cutoff, default = 0.05.
            
        Returns
        -------
        funpipe.plink
            An updated plink object with quality control file generated.

        """
        self._qc = self._bfile+'.qc'
        run(" ".join([
            'plink2 --bfile', self._bfile,
            '--recode --out '+ self._qc ,
            '--maf', maf,
            '--geno', miss,
            '--mind', ind
        ]))
        return self

#    def variant_qc(self):
        #         run("plink --bfile --test-missing --allow-extra-chr")
        # plink --bfile ${PREFIX} --cluster missing --allow-extra-chr
        # plink --bfile ${PREFIX} --missing --allow-extra-chr
#        return self

#    def sample_qc(self):
#        return self
