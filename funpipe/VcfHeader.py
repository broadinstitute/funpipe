from __future__ import division
import re
import sys
from scipy.stats import binom_test

# to do: merge with vcf.py


class VcfHeader:
    """
    Vcf Header Class
    """

    def __init__(self, vcf_file):
        """constructor of VcfHeader object
        
        Get samples, variant caller, sample columns, contig IDs and
        whether snpeff has been called.
        
        Parameters
        ----------
        vcf_file: string
            name of input vcf file
            
        
        """
        self.samples = []
        self.caller = None
        self.snpeff = False
        self.sample_columns = {}
        self.contigs = []

        comment_pattern = re.compile(r"^#")

        with open(vcf_file) as file:
            for full_line in file:
                line = full_line.rstrip()
                if re.search(comment_pattern, line):
                    if re.match('#CHROM', line):
                        fields = line.split('\t')
                        for i in range(9, len(fields)):
                            self.samples.append(fields[i])
                            self.sample_columns[fields[i]] = i
                    elif re.match('##PILON', line):
                        self.caller = 'PILON'
                    elif re.match('##GATK', line):
                        self.caller = 'GATK'
                    elif re.match('##SnpEff', line):
                        self.snpeff = True
                    elif re.match('##contig', line):
                        m = re.search('##contig=<ID=([^,]+),', line)
                        self.contigs.append(m.group(1))
                else:
                    break

        if self.samples == ['SAMPLE']:
            self.samples = [vcf_file]
            self.sample_columns[vcf_file] = 9
            sys.stderr.write("No sample name in " + vcf_file
                             + ", using file name.\n")

    def get_samples(self):
        return self.samples

    def get_caller(self):
        return self.caller

    def get_sample_index(self, sample):
        index = self.sample_columns[sample] - 9
        return index

    def get_snpeff_status(self):
        return self.snpeff

    def get_contigs(self):
        return self.contigs
