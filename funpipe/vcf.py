import os
from math import ceil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

import sys
#sys.path.append('.')
from funpipe.utils import run,rm
from funpipe.gt_pair import gt_pair
from funpipe.vcfheader import vcfheader
from funpipe.vcfrecord import vcfrecord

class vcf:
    def __init__(self, vcf_file, prefix='output', outdir='.', fasta=''):
        """Constructor of vcf object, that contains variant calling results.
        
        Parameters
        ----------
        vcf_file: string
            The path to input vcf file.
        prefix: string
            The output prefix, default = 'output'.
        outdir: string
            The output directory, default = '.'.
        fasta: string
            The path to input reference fasta file.
            
        Attributes
        ----------
        vcf: string
            The path to VCF file.
        vcfheader: funpipe.vcfheader
            The vcfheader object belonging to vcf object.
        prefix: string
            The prefix of output files.
        outdir: string
            The path to output directory.
        pairwise_share: np.ndarray()
            The matrix of # shared variants among all samples.
        pairwise_unique: np.ndarray()
            The matrix of # unique variants among all samples.
        n_samples: int
            The number of samples.
        dosage matrix: pd.DataFrame
            The genotype dosage matrix.
        af: pd.DataFrame
            The table of allelic frequency.
        ann_vcf: string
            The path to annotated VCF file.
        site_info: pd.DataFrame
            The table containing site level information.
        sample_info: pd.DataFrame
            The table containing sample level information.
        fasta: string
            The path to reference fasta file.
        
        """
        self._vcf = vcf_file
        self._vcfheader = vcfheader(self._vcf)
        #self._record_list = list()
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
        self._site_info_tsv = os.path.join( outdir, prefix + '.site_info.tsv')
        self._sample_info = pd.DataFrame()
        self._sample_info_tsv = os.path.join( outdir, prefix + '.sample_info.tsv' )
        self._fasta = fasta
        

    @property
    def site_info(self):
        """ Get site info dataframe.
        
        Returns
        -------
        pd.DataFrame
            The table containing site level information.
        """
        if self._site_info.empty:
            warnings.warn(
                "VCF site info dataframe is empty, get site info using 'get_site_info' function")
        return self._site_info
    
    @property
    def sample_info(self):
        """ Get sample info dataframe.
        
        Returns
        -------
        pd.DataFrame
            The table containing sample level information.
            
        """
        if self._sample_info.empty:
            warnings.warn(
                "VCF sample info dataframe is empty, get info using 'get_sample_info' function")
        return self._sample_info
    
#     @property
#     def record_list(self):
#         """Get the list of VCF records"""
#         if len(self._record_list) < 1:
#             warnings.warn(
#                     "The list of VCF records is empty, get VCF records using 'get_record' function")
#         return self._record_list
    
    @property
    def vcfheader(self):
        return self._vcfheader
    
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
        return self._var_counts

    @property
    def dosage_matrix(self):
        return self._dosage_matrix

    @staticmethod
    def create_snpeff_db(gff3, dest_dir, genome, config, prefix, ram=4, jar='/opt/snpEff/snpEff.jar', ref_fa = ''):
        """ Create snpEff database.
        
        Parameters
        ----------
        gff3: string
            The gff file of gene annotation.
        dest_dir: string
            The destination for snp database.
        genome: string
            The name of the reference genome.
        config: string
            The snpEff config files.
        prefix: string
            The output prefix.
        ram: int
            RAM in GB.
        jar: string
            The path to snpEff.jar.
        ref_fa: string
            The path to reference fasta file.
            
        Returns
        -------
        string
            Command line.
            
        """
        cmd = ' '.join(['create_snpeff_db.sh', dest_dir, jar, genome, ref_fa, gff3, str(ram) ])
        run(cmd)
        print("- SnpEff database created.")
        return cmd

    def snpeff_annot(self, jar='/opt/snpEff/snpEff.jar', config=None, genome=None, ram=None):
        """ Run SNPEFF on a VCF file.
        
        Parameters
        ----------
        jar: string
            The path to snpEff.jar
        config: string
            The snpEff config files.
        genome: string
            The tag of genome name.
        ram: int
            RAM in GB.
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with snpEff annotation generated.
            
        """
        self.ann_vcf = os.path.basename(self._vcf).replace('vcf', 'snpeff.vcf')
        run(' '.join([
            'java -Xmx'+str(ram)+'g', '-jar', jar, 'ann', '-v',
            '-c', config, '-i vcf -o vcf', genome,
            self._vcf, '| bgzip >', self.ann_vcf]))
        
        return self

    def import_snpeff(self, snpeff_tsv=None):
        """Import or generate snp site information, and organize into DataFrame.
        
        Parameters
        ----------
        snpeff_tsv: string
            Imported snp site information file, default = None.
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with a table of snp site level information generated.
            
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
        
        
        return self

    def af(self):
        """ Get allele frequencies using vcftools.
        
        Returns
        -------
        funpipe.vcf
            An updated vcf object with a table of allele frequencies generated.
            
        """
        run("vcftools --gzvcf "+self._vcf + " --freq2 --out tmp")
        header = ['chr','pos','n_alleles','n_chr','AF1','AF2']
        self._af = pd.read_csv('tmp.frq', sep='\t', header=None, names = header)
        self._af = self._af.drop(index = 0 )
        self._af = self._af.reset_index()
        self._af = self._af.drop(['index'],axis = 1)
        rm('tmp.frq')
        
        return self

    def get_sample_index(self, sample_id ):
        """ Get sample index from the VCF file.
        
        Parameters
        ----------
        sample_id: string
            The sample with index that the user looks for.
            
        Returns
        -------
        int
            The index of the sample.
            
        """
        header = vcfheader(self._vcf)
        sample_list = header.get_samples()
        
        if sample_list == [self._vcf]:
            return 0
        else:
            ind = sample_list.index(sample_id)
            
        return ind

    def select(self, jar='/opt/GATK-3.8/GenomeAnalysisTK.jar', outvcf=None, ref=None, snp=False, pass_only=False,
               indel=False):
        """ Parse VCF to get only sites passed quality control.

        Parameters
        ----------
        jar: string
            The path to GATK.jar.
        outvcf: string
            The path to output vcf.
        ref: string
            The path to reference.
        snp: bool
            Whether select snp variants.
        pass_only: bool
            If True, sites that are filtered will all be removed.
        indel: bool 
            Whether to select indel variants.

        Returns
        -------
        funpipe.vcf
            An updated vcf object with variants selected from the original VCF file.
            
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
        self._vcf = outvcf
        
        return self

    def filter_gt(self, outvcf, min_GQ=50, AD=0.8, DP=10):
        """ Apply variant filtering using GQ, AD and DP cutoffs.
        
        Parameters
        ----------
        outvcf: string
            The path to output filtered vcf.
        min_GQ: int
            Minimum GQ cutoff.
        AD: float
            Allelic depth cutoff.
        DP: int
            Depth cutoff.
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with genotypes filtered based on GQ,AD,DP in the original vcf file.
            
        """
        cmd = ' '.join(['filter_gatk_genotypes.py', '--min_GQ', str(min_GQ),
                        '--min_percent_alt_in_AD', str(AD),
                        '--min_total_DP', str(DP), self._vcf, '>', outvcf])
        self._vcf = outvcf
        run(cmd)
        
        return self

    def cal_dos(self, haploid=True):
        """ Get a genotype dosage matrix from the VCF file.
        
        Parameters
        ----------
        haploid: bool
            True if haploid, False if polyploid.
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with genotype dosage matrix generated.
            
        """
        dos_file = os.path.join( self._outdir, self._prefix+'.dos.tsv')
        if haploid:
            run("bcftools query -f '[%GT ]\\n' " + self._vcf + '>' + dos_file)
            self._dosage_matrix = pd.read_csv(
                dos_file, sep=r'\s+', header=None, na_values='.')
        else:
            raise ValueError("Not yet support polyploid.")
            
        return self

    def samples_concord(self, s1_idx, s2_idx, na_ignore=False):
        """ For each callset, compare SNP shared and unique between sample pairs.
        
        Parameters
        ----------
        s1_idx: string
            Sample 1's index in genotype dosages matrix.
        s2_idx: string
            Sample 2's index in genotype dosage matrix.
        na_ignore: bool
            Whether to ignore NaN data, default = False.
            
        Returns
        -------
        funpipe.gt_pair
            A genotype pair with computetd shared and unique variants.
            
        """
        gt = gt_pair(self._dosage_matrix[s1_idx], self._dosage_matrix[s2_idx],
                     na_ignore)
        
        gt.get_n_share()
        gt.get_n_unique()
        
        return gt

    def pairwise_concord(self, na_ignore=False):
        """ Get pairwise concordance amongst all sample pairs in the callset.
        
        Parameters
        ----------
        na_ignore: bool
            whether to ignore NaN data, default = False.
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with pairwise concordance computetd:
                * pairwise_share: numpy.ndarray, pairwise shared variants amongst all sample pairs.
                * pairwise_unique: numpy.ndarray, unique variants amongst all sample pairs.
                
        """
        if self._dosage_matrix is None:
            self.cal_dos()
            
        self._pairwise_share = np.zeros((self.n_samples, self.n_samples))
        self._pairwise_unique = np.zeros((self.n_samples, self.n_samples))
        for i in range(self.n_samples):
            for j in range(i, self.n_samples):
                gt = gt_pair(self._dosage_matrix[i], self._dosage_matrix[j],
                             na_ignore)
                
                self._pairwise_share[i][j] = gt.get_n_share()
                self._pairwise_unique[i][j] = gt.get_n_unique()
                
        return self

    
    def has_info(self, info, level = 'site' ):
        """Check whether an info field is present.
        
        Parameters
        ----------
        info: string
            Site or sample information.
        level: string
            The level = 'site' or 'sample', else raise value error, default = 'site'.
            
        Returns
        -------
        bool
            True, if the info field is present, else, False.
            
        """
        has_info = False
        sample_list = self._vcfheader.get_samples()
        
        if level == 'site':
            has_info = ( info in self._site_info.columns)
        elif level == 'sample':
            temp = True
            for i in range(len(sample_list)):
                temp = temp&( info+'_'+sample_list[i] in self._sample_info.columns)
            has_info = temp
        else:
            raise ValueError("Error, info level can only be site or sample.")
        
        if not has_info:
            raise ValueError( info+" is not presented in the "+level+" info column names." )
         
        return has_info
            
            
    # site level APIs
    def get_site_info(self, info=['AF']):
        """ Get variant site level infomation of interest.
        
        Parameters
        ----------
        info: list of string
            A list of site level information types, default = ['AF'].
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with a table of site level information generated:
                * site_info: pandas.dataframe, a table of site level information.
            
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
        return self
    
    
    def cal_maf(self, af_name='AF'):
        """Calculate minor allele frequency.
        
        Parameters
        ----------
        af_name: string
            The column name for allele frequencies
        
        Returns
        -------
        funpipe.vcf
            An updated vcf object with MAF computed and updated in site_info table.
            
        """
        self.has_info(af_name, 'site')
        self._site_info['MAF'] = self._af[af_name]
        self._site_info.ix[self._af[af_name] > 0.5, 'MAF'] = (
            1 - self._af.ix[self._af[af_name] > 0.5, 'MAF'])
        
        return  self

    def cal_miss(self):
        """ Calculate missingness of all sites.
        
        Returns
        -------
        funpipe.vcf
            An updated vcf obeject with missingness computed for all sites.
        
        """
        self.get_plink()
        run('plink2 --bfile '+self._plink+' --missing --allow-extra-chr --out '
            + os.path.join(self._outdir,self._prefix) )
        
        return self

    def plot_site_info(self, info, bins=100, density=True):
        """ Plot info field in the site. For example,minor allele frequencies.
        
        Parameters
        ----------
        info: string
            Info field in the site.
        bins: int
            The number of bins for histogram, default = 100.
        density: boolean
            If True, draw and return a probability density: each bin 
            will display the bin's raw count divided by the total number of counts.
            
        """
        self.has_info(info, 'site')

        plt.hist( self._site_info[info], bins=bins, density=density )
        plt.legend()
        plt.xlabel(info)
        plt.ylabel('# sites')

        return self

#     def genome_dist(self, info, pdf=None, window_size=2000000, ymax=12000, build='hg38',
#                     centro_tsv='/xchip/gtex/xiaoli/resources/centromeres.txt.tz'):
#         """plot distribution of genomic elements across the genome
        
#         Parameters
#         ----------
#         info: string
#             info field in the site.
#         pdf: string
#             output pdf name.
#         window_size: int
#             sliding windong size, default = 2000000.
#         ymax: int
#             y axis maximum, default = 12000
#         build: string
#             human genome build
#         centro_tsv: string
#             tsv containing centromere positions
            
#         """
#         has_info(info, 'site')
#         df = self._site_info
# #         centr = pd.read_csv(centro_tsv, header=None, sep='\t')
#         plt.figure(1)
#         for i in range(len(chrs)):
#             plt.subplot(4, 6, i+1)
#             pos = df.ix[df.contig == chrs[i], 'pos'].astype(float)
#             plt.ylim(0, ymax)
#             plt.hist( pos, bins=int(ceil(max(pos)/window_size)))
#             plt.title(chrs[i])
# #             centrStart = min(centr.loc[centr[0] == chrs[i], 1])
# #             centrEnd = max(centr.loc[centr[0] == chrs[i], 2])
# #             plt.axvline(x=centrStart, color='r')
# #             plt.axvline(x=centrEnd, color='r')
#         if pdf is not None:
#             plt.savefig(pdf)
#         else:
#             plt.show()
            
#         return(1)

    # sample level APIs
    def get_sample_info(self, info=['GT']):
        """ Get sample level info of interest
        
        Parameters
        ----------
        info: list of string
            A list of sample level information, default = ['GT','DP'].
            
        Returns
        -------
        funpipe.vcf
            An updated vcf object with a table of sample level information generated:
                * sample_info: pandas.dataframe, a table of sample level information.
            
        """
        sample_list = self._vcfheader.get_samples()
        header = ['CHR', 'POS','REF','ALT'] +  [ info[i]+'_'+sample_list[j] 
                                       for j in range(len(sample_list)) for i in range(len(info))]
        
        query_string = '\'%CHROM\t%POS\t%REF\t%ALT[\t'
        query_string += '\t'.join([('%'+i) for i in info])
        query_string =  query_string + ']\n\''
        
        cmd = ' '.join([
            "bcftools query -f ", query_string, self._vcf, '>',
            self._sample_info_tsv])
        run(cmd)
        self._sample_info = pd.read_csv(self._sample_info_tsv, sep='\t',
                                      header=None, names=header)
        return self
    
    #TODO
#     def plot_sample_info(self, info, bins=100, density=True, label=None ):
#         """ plot info field in the sample. For example,minor allele frequencies.
        
#         Parameters
#         ----------
#         info: string
#             info field in the sample.
#         bins: int
#             the number of bins for histogram, default = 100.
#         density: boolean
#             If True, draw and return a probability density: each bin 
#             will display the bin's raw count divided by the total number of counts.
#         label: string
#             plot labels.
            
#         """
        
#         has_info(info, 'sample')
        
#         plt.hist( self._sample_info[info], bins=bins, density=density, label = label )
#         plt.legend()
#         plt.xlabel(info)
#         plt.ylabel('# samples')
        
#         return self

    # formating
    def get_plink(self):
        """ Get plink format files 
        
        Returns
        -------
        funpipe.vcf
            an updated vcf object with bfile(bed,fam,bim) generated.
            
        """
        out = os.path.join(self._outdir,self._prefix)
        cmd = ' '.join(['plink2 --vcf', self._vcf,'--max-alleles 2' ,'--allow-extra-chr', '--make-bed --out',out])
        run(cmd)
        self._plink = out
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
