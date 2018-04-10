#!/usr/bin/env python3
import re
import gzip
import sys
from pipeline import cd, samtools_index_fa

'''
    Patch contig names in GFF and fasta from NCBI
    author: Xiao Li (xiaoli@broadinstitute.org)
'''

def chr_map_from_desc(desc_line):
    ''' extract contig and chromosome name from description line '''
    if 'mitochondrion' in desc_line:
        is_mt = True
        chr = re.sub(r'mitochondrion', 'chrMT', chr)
    elif 'chromosome':
        is_mt = False
        matched = re.match(r'>(\S+) .+ chromosome (\d+),.+', desc_line)
        contig, chr = matched.group(1), matched.group(2)
        if not chr:
            raise ValueError("contig name not matched!")
        chr = 'chr'+chr
    else:
        raise ValueError("Unrecognized description line:" + desc_line)
    return (contig, chr)


def patch_fasta(in_fa, out_fa, rm_mt=True):
    '''
    in_fa: input fasta file, can be gziped with suffix gz
    out_fa: output fasta file
    rm_mt: boolean, whether to remove mitochondrion sequence or not
    '''
    chr_map = {}
    if in_fa[-2:] == 'gz':
        fa = gzip.open(in_fa, 'rt')
    else:
        fa = open(in_fa, 'r')
    with open(out_fa, 'w') as pfa:
        for line in fa:
            if line[0] == '>':
                (contig, chr) = chr_map_from_desc(line)
                chr_map[contig] = chr
                line = '>'+chr+'\n'
            if rm_mt and (chr != 'chrMT'):
                pfa.write(line)
            else:
                pfa.write(line)
    return chr_map


def patch_gff_contig(chr_map, in_gff, out_gff, rm_mt=True):
    '''
    chr_map: a dictionary that maps
    in_gff: input gff file
    out_gff: output gff file
    rm_mt: boolean, whether remove mitochondrion genes or not
    '''
    with open(in_gff, 'r') as gff, open(out_gff, 'w') as pgff:
        for line in gff:
            if line[0] == '#':
                if line[:5] == '##seq':
                    seq_region = line.split(' ')
                    seq_region[1] = chr_map[seq_region[1]]
                    if rm_mt:
                        if seq_region[1] != 'mitochondrion':
                            pgff.write(' '.join(seq_region))
                    else:
                        pgff.write(' '.join(seq_region))
                else:
                    pgff.write(line)
            else:
                gff_line = line.split('\t')
                gff_line[0] = chr_map[gff_line[0]]
                if rm_mt:
                    if gff_line[0] != 'chrMT':
                        pgff.write('\t'.join(gff_line))
                else:
                    pgff.write('\t'.join(gff_line))

def write_chr_map(chr_map, out_tsv):
    with open(out_tsv, 'w') as out:
        for (k, v) in chr_map.items():
            out.write(k+'\t'+v+'\n')


if __name__ == '__main__':
    cd('/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21')
    prefix = 'GCF_000091045.1_ASM9104v1_genomic'
    chr_map = patch_fasta(
        prefix+'.fna.gz', prefix+'.patched.noMT.fasta', rm_mt=True)
    samtools_index_fa(prefix+'.patched.noMT.fasta')
    write_chr_map(chr_map, prefix+'_chr_map.tsv')
    # patch contig name in the chr
    patch_gff_contig(chr_map, prefix+'.gff', prefix+'.patched_noMT.gff')
