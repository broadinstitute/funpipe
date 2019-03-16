""" fungal analysis with a VCF file """
import os
from math import ceil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from funpipe.utils import run
import subprocess as sp


def pilon(fa, bam, prefix, ram, threads, jar):
    """ Run pilon commands

    Parameters
    ----------
        fa: :obj:`str` fasta file
        bam: :obj:`str` input bam path
        prefix: :obj:`str` output prefix
        ram: :obj:`int` input ram
        threads: :obj:`int` threads for pilon
        outdir: :obj:`str` output directory

    Returns
    -------


    """
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        '--genome', fa,
        '--frags', bam,
        '--output', prefix,
        '--threads', str(threads),
        '--vcf --changes --tracks --verbose > '+prefix+'.pilon.log 2>&1'])
    run(cmd)
    return cmd


def process_pilon_out(log, outdir, prefix):
    """ process pilon output
        log: logfile
        outdir: output directory
    """
    cmd = ' '.join(
         ['pilon_metrics', '-d', outdir, '-l', log, '--out_prefix', prefix])
    run(cmd)
    return cmd


def snpeff(invcf, outvcf, jar, config, genome, ram):
    """ run SNPEFF on a vcf
    invcf: input vcf
    outvcf: output vcf
    jar: snpeff jar
    genome: tag of genome name
    ram: memory in GB
    config: configuration file
    """
    cmd = ' '.join([
        'java -Xmx'+str(ram)+'g',
        '-jar', jar,
        'eff', '-v',
        '-c', config,
        '-onlyCoding False',
        '-i vcf -o vcf', genome, invcf, '>', outvcf])
    run(cmd)
    return cmd


def snpeff_db(gff3, dir, genome, config, prefix, ram, jar, ref_fa):
    """ Create snpEff database
    gff3: gff file of gene annotation
    genome: name of the reference genome
    config: snpEff config files
    prefix: output Prefix
    ram: RAM in GB
    jar: snpEff jar
    ref_fa: reference fasta file
    """
    snpeff_dir = os.path.dirname(jar)
    cmd = ' '.join(['sh snpeff_db.sh', dir, snpeff_dir, genome, ref_fa, gff3,
                    ram])
    run(cmd)
    return cmd


def tabix(file, type=None):
    """ Index tabix file
    :param file: input file
    :param type: file type, vcf
    """
    cmd = 'tabix '+file
    if type:
        cmd += ' -p '+type
    run(cmd)
    return file+'.tbi'


def filterGatkGenotypes(vcf, out_prefix):
    """ filter Gatk output vcf
    :param vcf: input vcf file
    :param out_prefix: output prefix
    """
    outfile = out_prefix+'_GQ50_AD08_DP10.vcf'
    cmd = ' '.join([
        'filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8',
        '--min_total_DP 10', vcf, '>', outfile
    ])
    run(cmd)
    return outfile


def filter_variants(invcf, outvcf, min_GQ=50, AD=0.8, DP=10):
    """ apply variant filtering using GQ, AD and DP
    :param invcf: input vcf
    :param outvcf: output vcf
    :param min_GQ: minimum GQ cutoff
    :param AD: allelic depth cutoff
    :param DP: depth cutoff
    """
    cmd = ' '.join(['filterGatkGenotypes.py', '--min_GQ', str(min_GQ),
                    '--min_percent_alt_in_AD', str(AD),
                    '--min_total_DP', str(DP), invcf, '>', outvcf])
    run(cmd)
    return outvcf


class vcf:
    """ vcf class command """
    def __init__(self, vcf, prefix='output', outdir='.'):
        self._vcf = vcf
        self._prefix = prefix
        self._outdir = outdir
        self._pairwise_share = None
        self._pairwise_unique = None
        self._var_counts = None
        self.n_samples = int(sp.check_output(
            'bcftools query -l '+vcf+' | wc -l', shell=True).decode().strip())
        self._dosage_matrix = None

    @property
    def nsamples(self):
        return self.n_samples

    @property
    def pairwise_concord(self):
        return self._pairwise_concord

    @property
    def pairwise_unique(self):
        return self._pairwise_unique

    @property
    def var_counts(self):
        return self.var_counts

    @property
    def dosage_matrix(self):
        return self._dosage_matrix

    def get_sample_index():
        """ return sample index from the VCF """
        return 0

    def select(self, jar, outvcf, ref, snp=False, pass_only=False,
               indel=False):
        """ parse VCF to get only sites passed quality control

        Parameters
        ----------
        jar: str
            GATK jar path
        prefix: str
            output file prefix
        outvcf:
            output vcf

        """
        if snp and indel:
            raise ValueError("Cannot select both SNPs and InDels")
        if not any([snp, pass_only, indel]):
            raise ValueError("At least select one type of variants")

        cmd = ' '.join([
            'java -jar ', jar, '-T SelectVariants', '-R', ref, '-V', self._vcf,
            '-o', outvcf
        ])
        cmd += ' -ef' if pass_only else ''
        cmd += ' -selectType SNP' if snp else ''
        cmd += ' -selectType INDEL' if indel else ''
        run(cmd)
        return self

    def filter_gt(self, outvcf, min_GQ=50, AD=0.8, DP=10):
        """ apply variant filtering using GQ, AD and DP
        :param invcf: input vcf
        :param outvcf: output vcf
        :param min_GQ: minimum GQ cutoff
        :param AD: allelic depth cutoff
        :param DP: depth cutoff
        """
        cmd = ' '.join(['filter_gatk_genotypes.py', '--min_GQ', str(min_GQ),
                        '--min_percent_alt_in_AD', str(AD),
                        '--min_total_DP', str(DP), self._vcf, '>', outvcf])
        self._vcf = outvcf
        run(cmd)
        return self

    def cal_dos(self, haploid=True):
        """ Get a genotype dosage matrix from the VCF """
        dos_file = self._prefix+'.dos.tsv'
        if haploid:
            run("bcftools query -f '[%GT ]\\n' " + self._vcf + '>' + dos_file)
            self._dosage_matrix = pd.read_csv(
                dos_file, sep=r'\s+', header=None, na_values='.')
        else:
            raise ValueError("Not yet support polyploid.")
        return self

    def samples_concord(self, s1_idx, s2_idx, na_ignore=False):
        """ For each callset, compare SNP sharing between sample pairs """
        gt = gt_pair(self._dosage_matrix[s1_idx], self._dosage_matrix[s2_idx],
                     na_ignore).get_n_unique()
        return gt.n_share, gt.n_unique

    def pairwise_concord(self, na_ignore=False):
        """ pairwise concordance amongst all sample pairs in the call set """
        if self._dosage_matrix is None:
            self = self.cal_dos()
        self._pairwise_share = np.zeros((self.n_samples, self.n_samples))
        self._pairwise_unique = np.zeros((self.n_samples, self.n_samples))
        for i in range(self.n_samples):
            for j in range(i, self.n_samples):
                gt = gt_pair(self._dosage_matrix[i], self._dosage_matrix[j],
                             na_ignore).get_n_unique()
                self._pairwise_share[i, j] = gt.n_share
                self._pairwise_unique[i, j] = gt.n_unique
        return self


class siteinfo:
    """ A table contain site level information """
    def __init__(self, ):
        self._vcf = ''
        self._df = pd.DataFrame()
        self._tsv = ''
        return self

    @property
    def vcf(self):
        return self._vcf

    @property
    def df(self):
        return self._df

    @property
    def tsv(self):
        return self

    # to do
    def import_tsv(self, tsv):
        return 0

    def import_vcf(self, vcf, info=['AF', 'AN', 'AC']):
        """ Import info from a VCF

        Description
        -----------
        get vcf and AF and missingness from a VCF
        Caveat: This module assumes the VCF's coming from GATK, with AF as the
        field for allele frequencies, and AC for Allele Count, and AN for
        Allelic Number.

        Parameters
        ----------
        VCF: str
            input VCF file path

        info: list
            A list that contains names of infor field of interest

        """
        header = ['CHR', 'ID'] + info
        query_string = '\'%CHROM\t%CHROM-%POS-%REF-%ALT{0}\t'
        query_string += '\t'.join([('%'+i) for i in info])+'\''

        cmd = ' '.join([
            "bcftools query -f ", query_string, self._vcf, '>', self._tsv])
        run(cmd)
        self._df = pd.read_csv(out_tsv, sep='\t', header=None, names=header)
        return self

    # def maf(self, ):
    #     """ calculate minor allele frequences """
    #
    # def miss(self,):
    #     """ calculate missingness """
    # # to do
    # def export_tsv(self,):
    #     return

    def cal_maf(self, af_name='AF'):
        """ calculate MAF
        :param df: data.frame containing allele frequencies
        :param AFname: column name for allele frequencies
        :rtype pandas dataframe
        """
        self._df['MAF'] = self._df[af_name]
        self._df.ix[self._df[af_name] > 0.5, 'MAF'] = (
            1 - self._df.ix[self._df[af_name] > 0.5, 'MAF'])
        return self

    def export_tsv(self, tsv):
        self._df.to_csv(tsv, sep='\t', index=False)
        return self

    # def dist_contrast(vec1, vec2, xlabel, ylabel, labels, pdf, bins=100):
    #     """ contract two distributions
    #     :param vec1: vector 1
    #     :param vec2: vector 2
    #     :param xlabel: label of x-axis
    #     :param ylable: label of y-axis
    #     :param pdf: output pdf name
    #     :param bins: number of bins for histograms
    #     """
    #     plt.figure()
    #     hist, bins = np.histogram(vec1, bins=100)
    #     plt.bar(bins[:-1], hist.astype(np.float32)/hist.sum(),
    #             width=(bins[1]-bins[0]), color='blue', alpha=0.5, label=labels[0])
    #     hist, bins = np.histogram(vec2, bins=100)
    #     plt.bar(bins[:-1], hist.astype(np.float32)/hist.sum(),
    #             width=(bins[1]-bins[0]), color='green', alpha=0.5, label=labels[1])
    #     plt.legend()
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    #     plt.savefig(pdf)
    #     plt.close()
    #     return(1)
    #
    #
    # def split_var_id(df, id_column='Variant'):
    #     """ split variant id column in a pandas data.frame
    #     :param df: pandas dataframe containing variant IDs
    #     :param id_column: column name of variant IDs
    #     """
    #     df[['contig', 'pos', 'ref', 'alt']] = df[id_column].str.split(':',
    #                                                                   expand=True)
    #     df = df.drop(id_column, axis=1)
    #     return(df)
    #
    #
    # def genome_dist(df, pdf=None, window_size=2000000, ymax=12000, built='hg38',
    #                 centro_tsv='/xchip/gtex/xiaoli/resources/centromeres.txt.tz'):
    #     """ distribution of genomic elements across the genome
    #     :param df: pandas dataframe containing
    #     :param pdf: output pdf name
    #     :param centro_tsv: tsv containing centromere positions
    #     :param window_size: sliding windong size
    #     :param ymax: y axis maximum
    #     :param build: human genome build
    #     :rtype boolean
    #     """
    #     centr = pd.read_csv(centro_tsv, header=None, sep='\t')
    #     if built == 'hg38':
    #         chrs = ['chr{0}'.format(i) for i in range(1, 23)]
    #         chrs.append('chrX')
    #     else:
    #         chrs = [format(i) for i in range(1, 23)]
    #         chrs.append('X')
    #     plt.figure(1)
    #     for i in range(len(chrs)):
    #         plt.subplot(4, 6, i+1)
    #         pos = df.ix[df.contig == chrs[i], 'pos'].astype(float)
    #         plt.ylim(0, ymax)
    #         pos.hist(bins=int(ceil(max(pos)/window_size)))
    #         plt.title(chrs[i])
    #         centrStart = min(centr.loc[centr[0] == chrs[i], 1])
    #         centrEnd = max(centr.loc[centr[0] == chrs[i], 2])
    #         plt.axvline(x=centrStart, color='r')
    #         plt.axvline(x=centrEnd, color='r')
    #     if pdf is not None:
    #         plt.savefig(pdf)
    #     else:
    #         plt.show()
    #     return(1)
    #
    #
    # def maf_plot(df, label='plot'):
    #     """ plot minor allele frequencies
    #     :param df: pandas dataframe containing VCF minor allele frequence fields
    #     :param label: plot labels
    #     """
    #     plt.figure()
    #     plt.hist(df.MAF,  bins=100, normed=True, label=label)
    #     plt.legend()
    #     plt.xlabel('MAF')
    #     plt.ylabel('# sites')
    #     plt.savefig(label+'.pdf')
    #     plt.close()
    #     return(1)


def _gt_type(gt):
    gt_type = type(gt).__name__
    if gt_type == 'Series':
        return gt
    elif gt_type in ['list', 'ndarray']:
        return pd.Series(gt_type)
    else:
        raise ValueError("Input gt vector should be either pandas series, list or numpy ndarray.")


class gt_pair:
    def __init__(self, gt1, gt2, na_ignore=False):
        """
        Parameters
        ----------
        gt1, gt2: pd.Series

        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt = gt_pair(gt1, gt2).get_n_unique()
        >>> print(gt.n_total, gt.n_unique, gt.n_share)
        7 5 2
        >>> gt = gt_pair(gt1, gt2, na_ignore=True).get_n_unique()
        >>> print(gt.n_total, gt.n_unique, gt.n_share)
        5 3 2

        """
        self.gt1 = _gt_type(gt1)
        self.gt2 = _gt_type(gt2)
        self.na_ignore = na_ignore
        self.n_total = None
        self.n_share = None
        self.n_unique = None
        self._not_both_ref = None

    def get_n_total(self):
        """ total number of non-monomorphic sites between two samples

        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_total().n_total
        7
        >>> gt_pair(gt1, gt2).get_n_total().n_total
        5

        """
        if self.na_ignore:
            self.n_total = ((self.gt1+self.gt2).fillna(0)
                            .map(lambda x: 1 if x !=0 else 0).sum())
        else:
            self.n_total = ((self.gt1.fillna(0) + self.gt2.fillna(0))
                            .map(lambda x: 1 if x !=0 else 0).sum())
        return self

    def get_n_share(self):
        """ Compare genotypes between two columns within a VCF, and report shared
        variants between the two samples.

                           A B
        for example: site1 1 .
                     site2 . 1
                     site3 1 1

        The unique variants here will be 2 (site1 and site2).

        Returns
        -------
        int: # shared sites

        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_share().n_share
        2

        Note
        ----

        This method is also cross-validated with GenotypeConcordance in GATK and
        bcftools stats.

        NaN will not be matched to any others.

        """
        # is a polymorphic site (non-reference sites)
        is_poly = (self.gt1 + self.gt2).map(lambda x: 1 if x != 0 else 0)
        # two sites are similar, include reference
        is_same = (self.gt1 - self.gt2).map(lambda x: 1 if x == 0 else 0)

        # number of shared alleles
        self.n_share = int((is_poly * is_same).sum())
        return self


    def get_n_unique(self):
        """
        Unique variants here mean a site that are private to either sample.

                               A B
        for example: site1 1 .
                     site2 . 1
                     site3 1 1

        The unique variants here will be 2 (site1 and site2). If ignore NA,
        the unique variants will be 0 (site1 and 2 will not be considered here).

        Parameters
        ----------
        gt1, gt2: pd.Series
        na_ignore:
            bool whether ignore na in the comparison

        Returns
        -------
            int: # unique sites

        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_unique()
        5
        >>> gt_unique(gt1, gt2)
        3

        """
        if self.n_total is None:
            self = self.get_n_total()
        if self.n_share is None:
            self = self.get_n_share()
        self.n_unique = self.n_total - self.n_share
        return self
