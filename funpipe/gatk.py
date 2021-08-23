import sys
sys.path.append('.')
import os
from utils import run


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
        fa: string
            input Fasta
        jar: string
            input jar file
        prefix: string
            output file prefix, default = \'output\'
        out_dir: string
            output directory for gatk commands, default = \'.\'
        RAM: int
            maximum RAM usage in gigabytes, default = 4
            
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
        vcf: string
            vcf: input vcf
        titv: bool      
            use TiTv Evaluator
        samp: bool
            stratify by samples
        indel: bool  
            use InDel Evaluator
        multi: bool
            summarize multiallelic sites
            
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
        vcf_dict: dict
            dictionary of vcf files, with key abbreviation of each vcf
        option: string
            merging options, UNIQUIFY,PRIORITIZE or UNSORTED.
        priority: 
            ordered list specifying priority for merging, default = None
            
        Returns
        -------
        string
            output combined vcf file
            
        '''
        out_vcf = self.out_dir + '/' +self.prefix+'.vcf.gz'
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
        in_vcf: string
            input vcf
        xl: string
            intervals to exclude
        il: string
            intervals to include
            
        Returns
        -------
        string
            output vcf file after selection
            
        '''
        output = self.out_dir+'/'+self.prefix+'.vcf.gz'
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
        comp_vcf: string
            VCF file for comparison
        eval_vcf: string
            VCF file for evaluation
        hap: bool
            whether input is haploid VCF
            
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
