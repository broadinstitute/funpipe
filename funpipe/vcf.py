""" 
Fungal analysis with a VCF file
===============================
"""
import os
from math import ceil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from funpipe.utils import run
from funpipe.gt_pair import gt_pair
import subprocess as sp


class vcf:
    """ vcf class command """
    def __init__(self, vcf_file, prefix='output', outdir='.', fasta=''):
        """variant call format
        
        Parameters
        ----------
        vcf_file: string
            input vcf file
        prefix: string
            output prefix, default = 'output'
        outdir: string
            output directory, default = '.'
        fasta: string
            input fasta file
            
        """
        self._vcf = vcf_file
        self._prefix = prefix
        self._outdir = outdir
        self._pairwise_share = None
        self._pairwise_unique = None
        self._var_counts = None
        self._n_samples = int(
            sp.check_output(
                'bcftools query -l '+vcf_file+' | wc -l', shell=True)
            .decode().strip())
        self._dosage_matrix = None
        self._af = None
        self.ann_vcf = None
        self._site_info = pd.DataFrame()
        self._site_info_tsv = prefix + '.site_info.tsv'
        self._sample_info = pd.DataFrame()
        self._sample_info_tsv = prefix + '.sample_info.tsv'
        self._fasta = fasta

    @property
    def info(self):
        """ Get site info dataframe"""
        if self._site_info.empty():
            warnings.warn(
                "VCF info dataframe is empty, get info using 'get_info' function")
        return self._site_info

    @property
    def n_samples(self):
        return self._n_samples

    @property
    def pairwise_share(self):
        return self._pairwise_share

    @property
    def pairwise_unique(self):
        return self._pairwise_unique

    @property
    def var_counts(self):
        return self.var_counts

    @property
    def dosage_matrix(self):
        return self._dosage_matrix

    @staticmethod
    def create_snpeff_db(gff3, dir, genome, config, prefix, ram, jar, ref_fa):
        """ Create snpEff database
        
        Parameters
        ----------
        gff3: string
            gff file of gene annotation
        dir: string
            destination for snpdb
        genome: string
            name of the reference genome
        config: string
            snpEff config files
        prefix: string
            output Prefix
        ram: int
            RAM in GB
        jar: string
            snpEff jar
        ref_fa: string
            reference fasta file
            
        Returns
        -------
        string
            command line
            
        """
        cmd = ' '.join(['create_snpeff_db.sh', dir, jar, genome, ref_fa, gff3, str(ram) ])
        run(cmd)
        return cmd

    def snpeff_annot(self, jar, config, genome, ram):
        """ run SNPEFF on a vcf
        
        Parameters
        ----------
        jar: string
            path to snpeff jar
        config: string
            configuration file
        genome: string
            tag of genome name
        ram:  int
            RAM in GB
            
        Returns
        -------
        string
            output snpeff file, annotated vcf
            
        """
        self.ann_vcf = os.path.basename(self._vcf).replace('vcf', 'snpeff.vcf')
        run(' '.join([
            'java -Xmx'+str(ram)+'g', '-jar', jar, 'ann', '-v',
            '-c', config, '-i vcf -o vcf', genome,
            self._vcf, '| bgzip >', self.ann_vcf]))
        return self.ann_vcf

    def import_snpeff(self, snpeff_tsv=None):
        """import or generate snp site information,
        and organize into DataFrame.
        
        Parameters
        ----------
        snpeff_tsv: string
            imported snp site information file, default = None
            
        Returns
        -------
        dataframe
            table of snp site level information
            
        """
        if snpeff_tsv is None:
            info_fields = [
                'AF', 'AN', 'AC',
                'SNPEFF_AMINO_ACID_CHANGE',
                'SNPEFF_CODON_CHANGE',
                'SNPEFF_EFFECT',
                'SNPEFF_EXON_ID',
                'SNPEFF_FUNCTIONAL_CLASS',
                'SNPEFF_GENE_BIOTYPE',
                'SNPEFF_GENE_NAME',
                'SNPEFF_IMPACT',
                'SNPEFF_TRANSCRIPT_ID'
            ]
            snpeff_tsv = self._prefix+'.snpeff.tsv'
            query = ('\'%CHROM\t%POS\t%REF\t%ALT\t'
                     + '\t'.join(['%INFO/'+i for i in info_fields])
                     + '\n\'')
            run('bcftools query -f {} '.format(query)
                + self._vcf+'> '+snpeff_tsv)
            
        self._site_info = pd.read_csv(
            snpeff_tsv, sep='\t', header=None,
            names=['CHR', 'POS', 'REF', 'ALT']+info_fields)
        
        
        return self._site_info

    def af(self):
        """ get allele frequencies using vcftools
        
        Returns
        -------
        dataframe
            table of allele frequencies
            
        """
        run("vcftools --gzvcf "+self._vcf + " --freq2 --out tmp")
        self._af = pd.read_csv('tmp.frq', sep='\t', header=0)
        rm('tmp.frq')
        return self._af

    def get_sample_index():
        """ return sample index from the VCF """
        return 0

    def select(self, jar, outvcf, ref, snp=False, pass_only=False,
               indel=False):
        """ parse VCF to get only sites passed quality control

        Parameters
        ----------
        jar: string
            GATK jar path
        outvcf: string
            output vcf
        ref: string
            reference
        snp: bool
            whether select snp variants
        pass_only: bool
            sites that are filtered will all be removed if True.
        indel: bool 
            whether select indel variants

        Returns
        -------
        string
            output selected vcf
            
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
        
        return outvcf

    def filter_gt(self, outvcf, min_GQ=50, AD=0.8, DP=10):
        """ apply variant filtering using GQ, AD and DP
        
        Parameters
        ----------
        outvcf: string
            output filtered vcf
        min_GQ: int
            minimum GQ cutoff
        AD: float
            allelic depth cutoff
        DP: int
            depth cutoff
            
        Returns
        -------
        string
            output filtered vcf
            
        """
        cmd = ' '.join(['filter_gatk_genotypes.py', '--min_GQ', str(min_GQ),
                        '--min_percent_alt_in_AD', str(AD),
                        '--min_total_DP', str(DP), self._vcf, '>', outvcf])
        self._vcf = outvcf
        run(cmd)
        
        return self._vcf

    def cal_dos(self, haploid=True):
        """ Get a genotype dosage matrix from the VCF
        
        Parameters
        ----------
        haploid: bool
            True if haploid, False if polyploid.
            
        Returns
        -------
        pandas.dataframe
            genotype dosage matrix
            
        """
        dos_file = self._prefix+'.dos.tsv'
        if haploid:
            run("bcftools query -f '[%GT ]\\n' " + self._vcf + '>' + dos_file)
            self._dosage_matrix = pd.read_csv(
                dos_file, sep=r'\s+', header=None, na_values='.')
        else:
            raise ValueError("Not yet support polyploid.")
            
        return self._dosage_matrix

    def samples_concord(self, s1_idx, s2_idx, na_ignore=False):
        """ For each callset, compare SNP sharing between sample pairs
        
        Parameters
        ----------
        s1_idx: string
            sample 1's index in genotype dosages matrix.
        s2_idx: string
            sample 2's index in genotype dosage matrix.
        na_ignore: bool
            whether ignore nan data, default = False.
            
        Returns
        -------
        int
            shared variants between 2 samples
        int
            unique variants between 2 samples
            
        """
        gt = gt_pair(self._dosage_matrix[s1_idx], self._dosage_matrix[s2_idx],
                     na_ignore)
        
        return gt.get_n_share(), gt.get_n_unique()

    def pairwise_concord(self, na_ignore=False):
        """ pairwise concordance amongst all sample pairs in the call set
        
        Parameters
        ----------
        na_ignore: bool
            whether ignore nan data, default = False
            
        Returns
        -------
        numpy.ndarray
            pairwise shared variants amongst all sample pairs
        numpy.ndarray
            unique variants amongst all sample pairs
        
        """
        if self._dosage_matrix is None:
            self.cal_dos()
            
        self._pairwise_share = np.zeros((self.n_samples, self.n_samples))
        self._pairwise_unique = np.zeros((self.n_samples, self.n_samples))
        for i in range(self.n_samples):
            for j in range(i, self.n_samples):
                gt = gt_pair(self._dosage_matrix[i], self._dosage_matrix[j],
                             na_ignore)
                
                self._pairwise_share[i, j] = gt.get_n_share()
                self._pairwise_unique[i, j] = gt.get_n_unique()
                
        return self._pairwise_share, self._pairwise_unique

    # site level APIs
    def has_info(self, info, type=['site', 'sample']):
        """Whether an info field is presented in the object
        
        Parameters
        ----------
        info: string
            site or sample information.
        type: string
            type = 'site' or 'sample', else raise value error.
            
        """
        has_info = False
        if type == 'site':
            has_info = (info in self._site_info.columns)
        elif type == 'sample':
            has_info = (info in self._sample_info.columns)
        else:
            raise ValueError("No such type")
        if not has_info:
            raise ValueError(
                info+" is not presented in the site info column names.")

    def get_info(self, info=['AF']):
        """ Get variant site level info of interest
        
        Parameters
        ----------
        info: list of string
            a list of site level information, default = ['AF'].
            
        Returns
        -------
        pandas.dataframe
            table of site level information
            
        """
        header = ['CHR', 'ID'] + info
        query_string = '\'%CHROM\t%CHROM-%POS-%REF-%ALT{0}\t'
        query_string += '\t'.join([('%'+i) for i in info])+'\''
        cmd = ' '.join([
            "bcftools query -f ", query_string, self._vcf, '>',
            self._site_info_tsv])
        run(cmd)
        self._site_info = pd.read_csv(self._site_info_tsv, sep='\t',
                                      header=None, names=header)
        return self._site_info
    
    
    def cal_maf(self, af_name='AF'):
        #TO DO
        """ calculate MAF
        
        
        :param df: data.frame containing allele frequencies
        :param AFname: column name for allele frequencies
        :rtype pandas dataframe
        """
        self.has_info(af_name, 'site')
        self._site_info['MAF'] = self._df[af_name]
        self._site_info.ix[self._df[af_name] > 0.5, 'MAF'] = (
            1 - self._df.ix[self._df[af_name] > 0.5, 'MAF'])
        
        return  self._site_info

    def cal_miss(self, name='miss'):
        """ Calculate missingness of all sites"""
        # TO DO
        self.get_plink()
        run('plink --bfile '+self._plink+' --missing --allow-extra-chr --out '
            +self._prefix)
        return self

    def plot_info(self, info, bins=100, normed=True):
        """ plot minor allele frequencies
        :param info: info field in the site
        :param df: pandas dataframe containing VCF minor allele frequence fields
        :param label: plot labels
        """
        self.has_info(info, site)

        plt.hist(df[info], bins=bins, normed=normed, label=label)
        plt.legend()
        plt.xlabel(info)
        plt.ylabel('# sites')

        return self

    def genome_dist(self, info, pdf=None, window_size=2000000, ymax=12000, built='hg38',
                    centro_tsv='/xchip/gtex/xiaoli/resources/centromeres.txt.tz'):
        """ distribution of genomic elements across the genome
        :param df: pandas dataframe containing
        :param pdf: output pdf name
        :param centro_tsv: tsv containing centromere positions
        :param window_size: sliding windong size
        :param ymax: y axis maximum
        :param build: human genome build
        :rtype boolean
        """
        has_info(info, 'site')

#         centr = pd.read_csv(centro_tsv, header=None, sep='\t')
        plt.figure(1)
        for i in range(len(chrs)):
            plt.subplot(4, 6, i+1)
            pos = df.ix[df.contig == chrs[i], 'pos'].astype(float)
            plt.ylim(0, ymax)
            pos.hist(bins=int(ceil(max(pos)/window_size)))
            plt.title(chrs[i])
#             centrStart = min(centr.loc[centr[0] == chrs[i], 1])
#             centrEnd = max(centr.loc[centr[0] == chrs[i], 2])
#             plt.axvline(x=centrStart, color='r')
#             plt.axvline(x=centrEnd, color='r')
        if pdf is not None:
            plt.savefig(pdf)
        else:
            plt.show()
        return(1)

    # sample level APIs
    def plot_sample_info(self, info):
        has_info(info, 'sample')
        return self

    # formating
    def get_plink(self):
        """ Get plink format files """
        cmd = ' '.join('plink --vcf', self._vcf, '--allow-extra-chr', '--out',
        self._prefix)
        run(cmd)
        self._plink = self._prefix
        return self

# incomplete class object
class siteinfo:
    """ A table contain site level information """
    def __init__(self):
        self._df = pd.DataFrame()
        self._vcf = ''
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

    def import_vcf(self, info=['AF', 'AN', 'AC']):
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
    #


# legacy methods for backward compatibility
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
