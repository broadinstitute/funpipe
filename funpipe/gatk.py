import os
import pandas as pd
from . import utils.run


class eval(analysis):
    def __init__(self, eval, out_dir, prefix, stats=None):
        """
        :param eval: input eval file f
        :param outdir:
        :param outfile: output file
        :param stats: Tables and statistics
        """
        if stats is None:
            stats = {
                'CountVariants': [
                    'nVariantLoci', 'variantRatePerBp', 'nSNPs',
                    'nInsertions', 'nDeletions', 'nHets', 'nHomRef', 'nHomVar',
                    'hetHomRatio', 'insertionDeletionRatio'
                ],
                'TiTvVariantEvaluator': ['nTi', 'nTv', 'tiTvRatio'],
                'IndelSummary': ['SNP_to_indel_ratio']
            }

    def parse_variant_eval(self):
        """ parse variantEval file
        :rtype
        """
        with open(self.eval, 'r') as fh:
            data = fh.read()
        tabs = data.rstrip().split('\n\n')

        meta_df = pd.DataFrame()
        for i in tabs:
            df = pd.read_csv(io.StringIO(i), comment='#', sep='\s+',
                             index_col='Sample')
            tab_name = df.columns[0]
            if tab_name in self.stats.keys():
                df = df[self.stats[tab_name]]
                meta_df = pd.concat([meta_df, df], axis=1)
        meta_df.to_csv(os.path.join(self.out_dir, self.prefix+'.tsv'),
                       sep='\t', compression='gzip')
        return meta_df


class gatk(analysis):
    """ Run GATK commands """
    def __init__(
            self, jar, prefix, out_dir='.', RAM=4,
            ):
        """
        :param jar: input GATK jar file
        """
        analysis.__init__(self, input, prefix='output', outdir='.', fasta=None,
                          gff=None, RAM=4, threads=1)
        self.cmd = ' '.join(['java -Xmx'+str(self.RAM)+'g -jar', jar,
                             '-R', self.fasta
        ])
        self.prefix = prefix

    def variant_eval(self, vcf, titv=True, samp=True, indel=True, multi=True):
        """ VCF sample QC by different stratifications
        :param vcf: input vcf
        :param titv: use TiTv Evaluator
        :param indel: use InDel Evaluator
        :param multi: summarize multiallelic sites
        :param samp: stratify by samples
        """
        out = os.path.join(self.out_dir, self.prefix+'.eval')
        cmd = ' '.join([self.cmd, '-T VariantEval', '--eval', vcf,
                        '-o', out, '-noEV -noST -EV CountVariants'])
        if titv:
            cmd += ' -EV TiTvVariantEvaluator'
        if samp:
            cmd += ' -ST Sample'
        if indel:
            cmd += ' -EV IndelSummary'
        if multi:
            cmd += ' -EV MultiallelicSummary'
        run(cmd)
        return out

    def combine_var(self, vcf_dict, option, priority=None):
        """
        :param vcf_dict: dictionary of vcf files, with key abbreviation of
                         each vcf
        :param prefix: output prefix
        :param option: merging options
        :param priority:
        """
        out_vcf = self.prefix+'.vcf.gz'
        options = ['UNIQUIFY', 'PRIORITIZE', 'UNSORTED']
        if option not in options:
            raise ValueError('Merge option not valid.\n')
        if option == 'PRIORITIZE' and priority is None:
            raise ValueError('Need to specify priority.\n')
        if option == 'UNSORTED':
            option += ' --assumeIdenticalSamples'
        cmd = ' '.join([
            self.cmd, '-T CombineVariants', '-genotypeMergeOptions', option,
            '-O', out_vcf])
        for name, vcf in vcf_dict.items():
            if option == 'PRIORITIZE':
                cmd += ' --variant:'+name+' '+vcf
            else:
                cmd += ' --variant '+vcf
        run(cmd)
        return out_vcf

    def select_var(self, in_vcf, xl=None, il=None):
        """ select variants
        :param in_vcf: input vcf
        :param xl: intervals to exclude
        :param il: intervals to include
        """
        output = self.prefix+'.vcf.gz'
        cmd = ' '.join([
            self.cmd, '-T SelectVariants', '--variant', in_vcf,
            '-o ', output])
        if xl is not None:
            cmd += '-XL '+xl
        if il is not None:
            cmd += '-L '+il
        run(cmd)
        return output


    def select_snps():


    def genotype_concordance(self, comp, eval):
        """ comppare
        :param comp: VCF file for comparison
        :parma eval: VCF file for evaluation
        :param out: output evaluation results
        """
        out = self.out_dir+'/'+self.prefix+'.txt'
        cmd = ' '.join([
            self.cmd, '-T GenotypeConcordance', '--comp', comp, '--eval', eval,
            '--out', out
        ])
        run(cmd)
        return out
