""" command with plink """


class plink:
    def __init__(self, prefix):
        self._bfile = prefix
        self._bed = prefix + '.bed'
        self._fam = prefix + '.fam'
        self._bim = prefix + '.bim'
        self._related = None

    def relatedness(self):
        """ Calculate genetic relatedness matrix """
        self._related = self._bfile+'.related.tsv'
        run(" ".join(['gemma', '-bfile '+self._bfile,
                      '-gk -o '+self._related]))
        return self

    def gwas(self, lmm=4):
        """ Fit a LMM for a univariate model with GEMMA"""
        self._assoc = self._bfile + '.gemma.assoc.tsv'
        run(" ".join([
            'gemma -bfile '+self._bfile,
            '-k '+self._related,
            '-lmm '+str(lmm),
            '-o '+self._assoc]))
        return self

    def import_pheno():
        """ import phenotypes """

    def gwas_filter(self, ind=0.1, miss=0.1, maf=0.05):
        """ Filter genotype and sample level missingness and AF

        Parameters
        ----------
            ind: individual level missingness
            miss: site level missingness
            maf: minor allele frequency cutoff

        """
        out_bfile = self._bfile+'.qc'
        run(" ".join([
            'plink --bfile', self._bfile,
            '--recode --out '+out_bfile,
            '--maf', maf,
            '--geno', miss,
            '--mind', ind
        ]))
        return self

    def variant_qc(self):
        #         run("plink --bfile --test-missing --allow-extra-chr")
        # plink --bfile ${PREFIX} --cluster missing --allow-extra-chr
        # plink --bfile ${PREFIX} --missing --allow-extra-chr
        return self

    def sample_qc(self):
        return self
