import sys
import os
from funpipe.utils import run

class gatk:
    def __init__(
            self, fa, jar='/opt/GATK-3.8/GenomeAnalysisTK.jar', prefix='output', out_dir='.', RAM=4):
        '''

        Parameters
        ----------
        fa: string
            The path to input Fasta.
        jar: string
            The path to GATK.jar.
        prefix: string
            The output file prefix, default = \'output\'.
        out_dir: string
            The output directory for gatk commands, default = \'.\'.
        RAM: int
            Maximum RAM usage in gigabytes, default = 4.
            
        Attributes
        ----------
        prefix: string
            The output file prefix, default = \'output\'.
        concordance: string
            The path to genotype concordance output file.
        eval: string
            The path to variant evaluation output file.
        combined_vcf: string    
            The path to combined variants output file.
        selected_vcf: string
            The path to selected variants output file.
            
        '''
        self.out_dir = out_dir
        self.cmd = ' '.join([
            'java -Xmx'+str(RAM)+'g -jar', jar, '-R', fa
        ])
        self.prefix = prefix
        self.concordance = None
        self.eval = None
        self.combined_vcf = None
        self.selected_vcf = None
          
    
    
    def variant_eval(self, vcf, titv=True, samp=True, indel=True, multi=True):
        ''' VCF sample QC by different stratifications.
        
        Parameters
        ----------
        vcf: string
            The path to input VCF file.
        titv: bool      
            Use TiTv Evaluator.
        samp: bool
            Stratify by samples.
        indel: bool  
            Use InDel Evaluator.
        multi: bool
            Summarize multiallelic sites.
            
        Returns
        -------
        funpipe.gatk
            An updated gatk object with the evaluation file of VCF generated,"prefix.eval".
            
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
        
        self.eval = out
        
        return self

    def combine_var(self, vcf_dict, option, priority=None):
        '''Combine variants from multiple VCF files.
        
        Parameters
        ----------
        vcf_dict: dict
            The input dictionary of VCF files, with key abbreviations of each VCF.
        option: string
            Merging options, UNIQUIFY,PRIORITIZE or UNSORTED.
        priority: 
            An ordered list specifying priority for merging, default = None.
            
        Returns
        -------
        funpipe.gatk
            An updated gatk object with combined VCF file generated,"prefix.combined.vcf.gz".
            
        '''
        out_vcf = os.path.join(self.out_dir,self.prefix+'.combined.vcf.gz')
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
        
        self.combined_vcf = out_vcf
        
        return self

    def select_var(self, in_vcf, xl=None, il=None):
        ''' Select variants by customized intervals to exclude and include.
        
        Parameters
        ----------
        in_vcf: string
            The path to input VCF file.
        xl: string
            Intervals to exclude.
        il: string
            Intervals to include.
            
        Returns
        -------
        funpipe.gatk
            An updated gatk object with the VCF file containing selected variants generated,"prefix.selected.vcf.gz".
            
        '''
        output = os.path.join(self.out_dir,self.prefix+'.selected.vcf.gz')
        cmd = ' '.join([
            self.cmd, '-T SelectVariants', '--variant', in_vcf,
            '-o ', output])
        if xl is not None:
            cmd += '-XL '+xl
        if il is not None:
            cmd += '-L '+il
        run(cmd)
        
        self.selected_vcf = output
        
        return self
    
    
    
    def genotype_concordance(self, comp_vcf, eval_vcf, hap=False):
        ''' Compare genotypes in 2 VCF files.
        
        Parameters
        ----------
        comp_vcf: string
            The path to VCF file for comparison.
        eval_vcf: string
            The path to VCF file for evaluation.
        hap: bool
            Whether input VCF files are haploid.
            
        Returns
        -------
        funpipe.gatk
            An updated gatk object with concordance output file generated,"prefix.concord.txt".
            
        '''
        out = os.path.join(self.out_dir,self.prefix+'.concord.txt')
        cmd = ' '.join([
            self.cmd, '-T GenotypeConcordance', '--comp', comp_vcf, '--eval', eval_vcf,
            '--out', out])
        run(cmd)
        
        self.concordance = out
        
        return self
 
