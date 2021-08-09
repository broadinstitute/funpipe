from __future__ import division
import re
import sys
from scipy.stats import binom_test

# to do: merge with vcf.py

class vcfrecord:
    """
    Object for each record (each line)
    """

    def __init__(self, vcf_line):
        """parse in each vcf record (each line) as
        chrom, pos, id, ref, alt, qual, filter, info, format
        and genotypes.
        
        Parameters
        ----------
        vcf_line: string
            a line of vcf record
        
        """
        self.vcf_line = vcf_line.rstrip()

        fields = self.vcf_line.split('\t')
        self.chrom = fields[0]
        self.pos = fields[1]
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = fields[7]
        self.format = fields[8]
        
        # problematic code
        self.genotypes = []
        self.vcf_annot = False
        for i in range(9, len(fields)):
            if fields[i]:
                gt_exp = re.compile(r"^[\d\.]")
                if re.search(gt_exp, fields[i]):
                    self.genotypes.append(fields[i])
                else:
                    self.vcf_annot = fields[i]

    def is_passing(self, caller):
        """ check whether the record passes Quality check
        
        Parameters
        ----------
        caller: string
            variant caller, either GATK or PILON
        
        Returns
        -------
        bool
            True if the record passess quality check, else False
            
        """
        if self.filter == 'PASS':
            return True
        elif caller == 'GATK' and self.filter == '.':  # not right
            return True
        else:
            return False

    def get_variant_type(self, caller, genotype):
        """ get variant type in vcf record
        
        Parameters
        ----------
        caller: string
            variant caller, either GATK or PILON
        
        genotype: string 
            genotype of an individual sample.
            The vertical pipe | indicates that the genotype is phased, 
            and is used to indicate which chromosome the alleles are on. 
            If this is a slash / rather than a vertical pipe, 
            it means we don’t know which chromosome they are on.
            
        Returns
        -------
        string
            variant type in [ 'uncalled_ambiguous', 'structural',
                'inside_deletion', 'SNP', 'DELETION', 'INSERTION',
                'unknown' ]
        """
        split_alt = self.alt.split(',')
        if genotype in ['.', './.', '.|.']:
            return 'uncalled_ambiguous'
        elif genotype in ['0', '0/0', '0|0']:
            return False
        else:
            alt = split_alt[int(genotype[-1:]) - 1]
            inequality_pattern = re.compile(r"^<\S+>$")
            if alt == '.':
                return False
            elif caller == 'PILON' and re.search(inequality_pattern, alt):
                return 'structural'
            elif caller == 'GATK' and re.search(inequality_pattern, alt):
                return 'inside_deletion'
            elif len(alt) == 1 and len(self.ref) == 1:
                return 'SNP'
            elif len(alt) < len(self.ref):
                return 'DELETION'
            elif len(alt) >= len(self.ref):
                return 'INSERTION'
            else:
                return 'unknown'

    def get_variant_length(self, genotype):
        """ Get Variant length for a specific sample
        
        Parameters
        ----------
        genotype: string
            genotype of an individual sample
            The vertical pipe | indicates that the genotype is phased, 
            and is used to indicate which chromosome the alleles are on. 
            If this is a slash / rather than a vertical pipe, 
            it means we don’t know which chromosome they are on.
            
        Returns
        -------
        int
            length of a variant
            
        """
        if genotype in ['0', '0/0', '0|0', '.', '.|.', './.']:
            return False
        else:
            ref_length = len(self.get_ref())
            genotype_end = re.search("(\d+)$", genotype)
            alt_length = len(self.get_alt(genotype_end.group(1)))
            if ref_length > 1 or alt_length > 1:
                return abs(alt_length - ref_length)
            else:
                return int(0)

    def get_genotype(self, index=0, min_gq=0, min_per_ad=float(0),
                     min_tot_dp=0, het_binom_p=False, return_flags=False):
        """Return genotypes of a record
        
        Parameters
        ----------
        index: int
            index of genotype, default = 0
        min_gq: int
            minimum genotype quality score cutoff, default = 0
        min_per_ad: float
            minimum allelic depth, default = float(0)
        min_tot_dp: int
            minimum total depth, default = 0
        het_binom_p: bool
            whether compute pvalue of heterozygous binomial test
        return_flags: bool
            if True, return parsed genotype, else
            return parsed genotype list.
            
        :return:
        """
        # working on function to accommodate hets
        genotype = self.genotypes[index]
        parsed_genotype = genotype.split(':')[0]
        dip_flag = False
        if len(parsed_genotype) == 3:
            dip_flag = True
        parsed_genotype_list = [parsed_genotype, 0, 0, 0, 0]
        # print(parsed_genotype) ####
        try:
            gq = self.get_GQ(parsed_genotype, index)
            if int(gq) < int(min_gq):
                parsed_genotype = '.'
                if dip_flag:
                    parsed_genotype = './.'
                parsed_genotype_list[1] = 1
        except:
            pass
        try:
            percent_ad = self.get_percent_AD(index)
            if float(percent_ad) < min_per_ad:
                parsed_genotype = '.'
                if dip_flag:
                    parsed_genotype = './.'
                parsed_genotype_list[2] = 1
        except:
            pass
        try:
            total_dp = self.get_total_DP(index)
            if int(total_dp) < int(min_tot_dp):
                parsed_genotype = '.'
                if dip_flag:
                    parsed_genotype = './.'
                parsed_genotype_list[3] = 1
        except:
            pass
        het_flag = self.is_het(index)
        if het_binom_p and het_flag:
            try:
                p = self.get_AD_binomial_p(index)
                if p < float(het_binom_p):
                    parsed_genotype = './.'
                    parsed_genotype_list[4] = 1
            except:
                pass
        if not return_flags:
            return parsed_genotype
        else:
            parsed_genotype_list[0] = parsed_genotype
            return parsed_genotype_list

    def is_het(self, index=0):  # currently works only on biallelic sites
        """check whether a genotype is heterogenous (currently works only on biallelic sites)
        
        Parameters
        ----------
        index: int
            index for GT field in genotype string
            
        Returns
        -------
        bool
            True if the genotype is heterogenous, else False.
        """
        het = False
        genotype = self.genotypes[index]
        parsed_genotype = genotype.split(':')[0]
        gt_list = list(parsed_genotype)
        if '1' in gt_list and '0' in gt_list:
            het = True
        return het

    def get_GQ(self, parsed_genotype, index=0):
        """get GQ(genotype quality)
        
        Parameters
        ----------
        parsed_genotype: string
            genotype parsed in
        index: int
            index of genotype, default = 0
            
        Returns
        --------
        string
            genotype quality
            
        """
        gq = 'Undefined'
        fields = self.genotypes[index]
        gt_fields = fields.split(':')
        try:
            gq_index = self.get_GQ_index(parsed_genotype)
            gq = gt_fields[gq_index]
        except:
            pass
        return gq

    def get_percent_AD(self, index=0):
        # currently works only on biallelic sites
        
        """get percentage of allelic depths,
        alt_depth/(alt_depth + ref_depth) 
        
        Parameters
        ----------
        index: int
            index of genotype, default = 0
        
        Returns
        -------
        float
            percentage of allelic depths
        
        Example
        -------
        For a GT field,  "GT:AD    0/0:6,9", percentage is 9/(6+9).
        
        """
        percent_AD = 'Undefined'
        fields = self.genotypes[index]
        gt_fields = fields.split(':')
        try:
            ad_index = self.get_AD_index()
            ad = gt_fields[ad_index]
            split_ad = ad.split(',')
            if split_ad[0] != '0':
                percent_AD = (float(split_ad[1]) / (float(split_ad[0])
                                                    + float(split_ad[1])))
            elif split_ad[1] != '0':
                percent_AD = float(1)
        except:
            pass
        return percent_AD

    def get_AD_binomial_p(self, index=0):
        # currently works only on biallelic sites
        
        """binomial test for allelic depths(ref and alt),
        assess the significance with a p-value of the hypothesis test.
        The hypothesized probability of success is 0.5, and the test is 2 sided.
        
        Parameters
        ----------
        index: int
            index of genotype, default = 0
            
        Returns
        -------
        float
            p value of binomial test
        
        """
        pvalue = None
        fields = self.genotypes[index]
        gt_fields = fields.split(':')
        try:
            ad_index = self.get_AD_index()
            ad = gt_fields[ad_index]
            split_ad = ad.split(',')
            if split_ad[1] != '0':
                pvalue = binom_test(int(split_ad[1]),
                                    (int(split_ad[0]) + int(split_ad[1])))
            elif split_ad[0] != '0':
                pvalue = binom_test(int(split_ad[0]),
                                    (int(split_ad[1]) + int(split_ad[0])))
        except:
            pass
        
        return pvalue

    def get_total_DP(self, index=0):
        """get total read depth(filtered)
        
        Parameters
        ----------
        index: int
            index of genotype, default = 0
            
        Returns
        -------
        string
            total read depth
            
        """
        total_DP = 'Undefined'
        fields = self.genotypes[index]
        gt_fields = fields.split(':')
        try:
            dp_index = self.get_DP_index()
            total_DP = gt_fields[dp_index]
        except:
            pass
        return total_DP

    def get_GQ_index(self, parsed_genotype):
        """get index of GQ() in format field
        
        Parameters
        ----------
        parsed_genotype: string
            genotype parsed in
            
        Returns
        -------
        int
            index of GQ
            
        Example
        -------
            AD:GQ:DP:HQ
        index of GQ is 1 in format field.
        
        """
        gq_index = 'Undefined'
        fields = self.format
        format_fields = fields.split(':')
        if parsed_genotype == '0' or parsed_genotype == '0/0' or parsed_genotype == '0|0':
            try:
                gq_index = format_fields.index('RGQ')
            except:
                try:
                    gq_index = format_fields.index('GQ')
                except:
                    pass
        else:
            try:
                gq_index = format_fields.index('GQ')
            except:
                pass
        return gq_index

    def get_AD_index(self):
        """get the index of AD in format field,
        AD is allelic depths for the ref and alt alleles in the order listed.
        
        Returns
        -------
        int
            index of AD
            
        Example
        -------
            AD:GQ:DP:HQ
        index of AD is 0 in format field.
        
        """
        ad_index = 'Undefined'
        fields = self.format
        format_fields = fields.split(':')
        try:
            ad_index = format_fields.index('AD')
        except:
            pass
        return ad_index

    def get_DP_index(self):
        """get the index of DP(valid read depth) in format field
        
        Returns
        -------
        int
            index of DP
            
        Example
        -------
            GT:GQ:DP:HQ
        index of DP is 2 in format field.
        
        """
        ad_index = 'Undefined'
        fields = self.format
        format_fields = fields.split(':')
        try:
            ad_index = format_fields.index('DP')
        except:
            pass
        return ad_index

    def get_chrom(self):
        return self.chrom

    def get_pos(self):
        return self.pos

    def get_id(self):
        return self.id

    def get_ref(self):
        return self.ref

    def get_alt(self, genotype):
        """get alternative allele
        
        Parameters
        ----------
        genotype: string
            input genotype in vcf record, for example, 0|0
        
        Returns
        -------
        string
            alternative allele
            
        """
        if genotype in ['0', '.', '0/0', './.', '0|0', '.|.']:
            return False
        else:
            split_genotype = re.split(r"[/|]", genotype)
            if len(set(split_genotype)) > 1:
                return 'N'
            else:
                split_alt = self.alt.split(',')
                genotype_end = re.search('(\d+)$', genotype)
                index = int(genotype_end.group(1)) - 1
                return split_alt[index]

    def get_alt_field(self):
        return self.alt

    def get_snpeff_annot(self, alt):
        fields = self.info.split(';')
        for field in fields:
            if re.match('ANN', field):
                ann_fields = field.split(',')
                # print ann_fields ##
                for ann_field in ann_fields:
                    m = re.match('(ANN=)*(\w+)', ann_field)
                    if m.group(2):
                        if m.group(2) == alt:
                            return ann_field.replace('ANN=', '')
        return False

    def get_snpeff_effect(self, snpeff_annot):
        try:
            annot_sections = snpeff_annot.split('|')
            return annot_sections[1]
        except:
            return False

    def get_snpeff_impact(self, snpeff_annot):
        try:
            annot_sections = snpeff_annot.split('|')
            return annot_sections[2]
        except:
            return False

    def get_snpeff_feature(self, snpeff_annot):
        try:
            annot_sections = snpeff_annot.split('|')
            return annot_sections[3]
        except:
            return False

    def get_qual(self):
        return self.qual

    def get_filter(self):
        return self.filter

    def get_info(self):
        return self.info

    def get_AF(self):
        """get allelic frequency at locus
        
        Returns
        -------
        float
            allelic frequency
        """
        AF = None
        fields = self.info.split(';')
        for field in fields:
            m = re.search('AF=([\d\.]+)', field)
            try:
                AF = float(m.group(1))
            except:
                raise Exception('Sorry, AF information not found in vcf record.')
                
        return AF

    def get_QP(self):
        """get QP: Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus
        
        Returns
        -------
        dict
            A,C,G,T weights at locus
        """
        QP_W = None
        fields = self.info.split(';')
        for field in fields:
            m = re.search('QP=(\d+),(\d+),(\d+),(\d+)', field)
            try:
                QP = m.group(1, 2, 3, 4)
                QP = [int(q) for q in QP]
                QP_W = {'A':QP[0],'C':QP[1],'G':QP[2],'T':QP[3]}
            except:
                raise Exception('Sorry, QP information not found in vcf record.')
                
        return QP_W

    def get_MAF_from_QP(self):
        MAF = False
        QP = self.get_QP()
        try:
            MAF = (100 - max(QP)) / 100
        except:
            pass
        return MAF

    def get_format(self):
        return self.format

    def get_genotypes_fields(self):
        return self.genotypes

    def get_genotypes_field(self, sample_index):
        return self.genotypes[sample_index]

    def get_vcf_annot(self):
        return self.vcf_annot

    def is_singleton(self):
        #ToDo
        """
        """
        nonzeros = 0
        last_nonzero_index = False
        for i in range(0, len(self.genotypes)):
            current_genotype = self.get_genotype(index=i)
            if re.search(r"[1-9]", current_genotype):
                nonzeros += 1
                last_nonzero_index = i
        if nonzeros == 1:
            return str(last_nonzero_index)
        else:
            return False

    def is_biallelic(self):
        """check whether the site is biallelic
        
        Returns
        -------
        bool
            True if the site is biallelic, else False
        """
        split_alt = self.alt.split(',')
        if len(split_alt) > 1:
            return False
        else:
            return True

    def count_ambig_genotypes(self):
        """count the number of ambiguous genotypes
        
        Returns
        -------
        int
            the number of ambiguous genotypes
        """
        ambig = 0
        for current_genotype in self.genotypes:
            if current_genotype.split(':')[0] in ['.', './.']:
                ambig += 1
        return ambig

    def get_genotype_profile(self):
        """get a profile of all samples' genotypes in this vcf record
        
        Returns
        -------
        list of strings
            genotype profile of all samples
        """
        profile = list()
        for current_genotype in self.genotypes:
            profile.append(current_genotype.split(':')[0])
        return profile


