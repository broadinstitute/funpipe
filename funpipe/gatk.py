import os
from funpipe.utils import run


"""
GATK: GenomeAnalysisToolkit
===========================
"""

class gatk:
    """ Run GATK commands """
    def __init__(
            self, fa, jar, prefix='output', out_dir='.', RAM=4):
        ''' constructor of gatk object
        
        Parameters
        ----------
        arg1: string
            fa: input Fasta
        arg2: string
            jar: input jar file
        arg3: string
            prefix: output file prefix, default = \'output\'
        arg4: string
            out_dir: output directory for gatk commands, default = \'.\'
        arg5: int
            RAM: RAM usage, default = 4
            
        '''
        self.out_dir = out_dir
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar, '-R', fa
        ])
        self.prefix = prefix

    def variant_eval(self, vcf, titv=True, samp=True, indel=True, multi=True):
        ''' VCF sample QC by different stratifications
        
        Parameters
        ----------
        arg1: string
            vcf: input vcf
        arg2: bool      
            titv: use TiTv Evaluator
        arg3: bool
            samp: stratify by samples
        arg4: bool  
            indel: use InDel Evaluator
        arg5: bool
            multi: summarize multiallelic sites
            
        Returns
        -------
        string
            output file name
            
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
        '''combine variants
        
        Parameters
        ----------
        arg1: dict
            vcf_dict: dictionary of vcf files, with key abbreviation of each vcf
        arg2: string
            option: merging options, UNIQUIFY,PRIORITIZE or UNSORTED.
        arg3: 
            priority:  default = None
            
        Returns
        -------
        string
            output combined vcf file
            
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
        
        Parameters
        ----------
        arg1: string
            in_vcf: input vcf
        arg2: string
            xl: intervals to exclude
        arg3: string
            il: intervals to include
            
        Returns
        -------
        string
            output vcf file after selection
            
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

    def genotype_concordance(self, comp_vcf, eval_vcf, hap=False):
        ''' comppare 2 vcf files
        
        Parameters
        ----------
        arg1: string
            comp_vcf: VCF file for comparison
        arg2: string
            eval_vcf: VCF file for evaluation
        arg3: bool
            hap: whether input is haploid VCF
            
        Returns
        -------
        string
            output comparison result name
            
        '''
        out = self.out_dir+'/'+self.prefix+'.txt'
        cmd = ' '.join([
            self.cmd, '-T GenotypeConcordance', '--comp', comp_vcf, '--eval', eval_vcf,
            '--out', out
        ])
        run(cmd)
        return out
