from math import ceil
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from . import utils.run


class vcf(analysis):
    def __init__(self, vcf, fa, prefix, outdir='.' RAM=4, threads=1):
        """ Initialize a VCF object

        Parameters
        ----------
            vcf: path to input vcf file
            fa: reference fasta file
            prefix: vcf file prefix
            outdir: output directory, default '.'
            RAM: RAM for the analysis
            threads: number of threads used for this analysis

        """
        analysis.__init__(self, input, prefix='output', outdir='.', fasta=None,
                          gff=None, RAM=4, threads=1)
        self.vcf = vcf
        self.fa = fa
        self.prefix = prefix
        self.outdir = outdir
        self.RAM = RAM
        self.threads = threads

    def process_pilon_out(self.log, outdir, prefix):
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

    def filter_genotypes(self, prefix, min_GQ=50, AD=0.8, DP=10):
        """ apply variant filtering using GQ, AD and DP
        :param prefix: output vcf prefix
        :param min_GQ: minimum GQ cutoff
        :param AD: allelic depth cutoff
        :param DP: depth cutoff
        """
        output_vcf = prefix+'.vcf'
        cmd = ' '.join(['filterGatkGenotypes.py',
                        '--min_GQ', str(min_GQ),
                        '--min_percent_alt_in_AD', str(AD),
                        '--min_total_DP', str(DP),
                        self.vcf, '>', output_vcf])
        run(cmd)
        return output_vcf

    def get_vcf_af_miss(vcf, out_tsv):
        """ get vcf and AF and missingness from a VCF
        Caveat: This module assumes the VCF's coming from GATK, with AF as the
        field for allele frequencies, and AC for Allele Count, and AN for Allele
        Number.
        :param vcf: input path of vcf
        :param out_tsv: output tsv file name
        """
        cmd = ''
        return out_tsv

    def cal_maf(df, af_name='AF'):
        """ calculate MAF
        :param df: data.frame containing allele frequencies
        :param AFname: column name for allele frequencies
        :rtype pandas dataframe
        """
        df['MAF'] = df[af_name]
        df.ix[df[af_name] > 0.5, 'MAF'] = 1 - df.ix[df[af_name] > 0.5, 'MAF']
        return df

    def dist_contrast(vec1, vec2, xlabel, ylabel, labels, pdf, bins=100):
        """ contract two distributions
        :param vec1: vector 1
        :param vec2: vector 2
        :param xlabel: label of x-axis
        :param ylable: label of y-axis
        :param pdf: output pdf name
        :param bins: number of bins for histograms
        """
        plt.figure()
        hist, bins = np.histogram(vec1, bins=100)
        plt.bar(bins[:-1], hist.astype(np.float32)/hist.sum(),
                width=(bins[1]-bins[0]), color='blue', alpha=0.5, label=labels[0])
        hist, bins = np.histogram(vec2, bins=100)
        plt.bar(bins[:-1], hist.astype(np.float32)/hist.sum(),
                width=(bins[1]-bins[0]), color='green', alpha=0.5, label=labels[1])
        plt.legend()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(pdf)
        plt.close()
        return(1)

    def split_var_id(df, id_column='Variant'):
        """ split variant id column in a pandas data.frame
        :param df: pandas dataframe containing variant IDs
        :param id_column: column name of variant IDs
        """
        df[['contig', 'pos', 'ref', 'alt']] = df[id_column].str.split(':',
                                                                      expand=True)
        df = df.drop(id_column, axis=1)
        return(df)

    def genome_dist(df, pdf=None, window_size=2000000, ymax=12000, built='hg38',
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
        centr = pd.read_csv(centro_tsv, header=None, sep='\t')
        if built == 'hg38':
            chrs = ['chr{0}'.format(i) for i in range(1, 23)]
            chrs.append('chrX')
        else:
            chrs = [format(i) for i in range(1, 23)]
            chrs.append('X')
        plt.figure(1)
        for i in range(len(chrs)):
            plt.subplot(4, 6, i+1)
            pos = df.ix[df.contig == chrs[i], 'pos'].astype(float)
            plt.ylim(0, ymax)
            pos.hist(bins=int(ceil(max(pos)/window_size)))
            plt.title(chrs[i])
            centrStart = min(centr.loc[centr[0] == chrs[i], 1])
            centrEnd = max(centr.loc[centr[0] == chrs[i], 2])
            plt.axvline(x=centrStart, color='r')
            plt.axvline(x=centrEnd, color='r')
        if pdf is not None:
            plt.savefig(pdf)
        else:
            plt.show()
        return(1)

    def maf_plot(df, label='plot'):
        """ plot minor allele frequencies
        :param df: pandas dataframe containing
        :param label: plot labels
        """
        plt.figure()
        plt.hist(df.MAF,  bins=100, normed=True, label=label)
        plt.legend()
        plt.xlabel('MAF')
        plt.ylabel('# sites')
        plt.savefig(label+'.pdf')
        plt.close()
        return(1)
