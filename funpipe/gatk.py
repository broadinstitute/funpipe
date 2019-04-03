import os
from funpipe.utils import run


class gatk:
    """ Run GATK commands """
    def __init__(
            self, fa, jar, prefix='output', out_dir='.', RAM=4):
        ''' VCF sample QC
        :param fa: input Fasta
        :param prefix: output file prefix
        :param jar: input jar file
        :param RAM: RAM usage
        :param out_dir: output directory
        '''
        self.out_dir = out_dir
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar, '-R', fa
        ])
        self.prefix = prefix

    def variant_eval(self, vcf, titv=True, samp=True, indel=True, multi=True):
        ''' VCF sample QC by different stratifications
        :param vcf: input vcf
        :param titv: use TiTv Evaluator
        :param indel: use InDel Evaluator
        :param multi: summarize multiallelic sites
        :param samp: stratify by samples
        '''
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
        '''
        :param vcf_dict: dictionary of vcf files, with key abbreviation of
                         each vcf
        :param prefix: output prefix
        :param option: merging options
        :param priority:
        '''
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
        ''' select variants
        :param in_vcf: input vcf
        :param xl: intervals to exclude
        :param il: intervals to include
        '''
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

    def genotype_concordance(self, comp, eval, hap=False):
        ''' comppare
        :param comp: VCF file for comparison
        :parma eval: VCF file for evaluation
        :param out: output evaluation results
        :param hap: whether input is haploid VCF
        '''
        out = self.out_dir+'/'+self.prefix+'.txt'
        cmd = ' '.join([
            self.cmd, '-T GenotypeConcordance', '--comp', comp, '--eval', eval,
            '--out', out
        ])
        run(cmd)
        return out
