#!/usr/bin/env python3
import re
import gzip
from pipeline import cd

'''
    Patch contig names in GFF and fasta from NCBI
    author: Xiao Li (xiaoli@broadinstitute.org)
'''

def patch_fasta(in_fa, out_fa, rm_mt = True):
    '''
    in_fa: input fasta file, can be gziped with suffix gz
    out_fa: output fasta file
    rm_mt: boolean, whether to remove mitochondrion sequence or not
    '''
    is_mt = False
    chr_map = {}
    if in_fa[-2:] == 'gz':
        fa = gzip.open(in_fa, 'rt')
    else:
        fa = open(in_fa, 'r')
    with open(out_fa, 'w') as pfa:
        for line in fa:
            if line[0] == '>':
                matched = re.match(
                    r'>(\S+) Cryptococcus neoformans var. grubii H99 (.+),.+',
                    line)
                contig, chr = matched.group(1), matched.group(2)
                if 'mitochondrion' in line:
                    is_mt = True
                    chr = re.sub(r'mitochondrion', 'chrMT', chr)
                else:
                    is_mt = False
                    chr = re.sub(r'chromosome ', 'chr', chr)
                chr_map[contig] = chr
                line = '>'+chr+'\n'
            if rm_mt and (not is_mt):
                pfa.write(line)
    return chr_map


def patch_gff_contig(chr_map, in_gff, out_gff):
    with open(in_gff, 'r') as gff, open(out_gff, 'w') as pgff:
        for line in gff:
            if line[0] == '#':
                if line[:5] == '##seq':
                    seq_region = line.split(' ')
                    seq_region[1] = chr_map[seq_region[1]]
                    pgff.write(' '.join(seq_region))
                else:
                    pgff.write(line)
            else:
                gff_line = line.split('\t')
                gff_line[0] = chr_map[gff_line[0]]
                pgff.write('\t'.join(gff_line))

def write_chr_map(chr_map, out_tsv):
    with open(out_tsv, 'w') as out:
        for (k, v) in chr_map.items():
            out.write(k+'\t'+v+'\n')


if __name__ == '__main__':
    cd('/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18')
    chr_map = patch_fasta(
        'GCF_000149245.1_CNA3_genomic.fna',
        'GCF_000149245.1_CNA3_genomic.patched.noMT.fna',
        rm_mt=True)
    write_chr_map(chr_map, 'GCF_000149245.1_CNA3_chr_map.tsv')
    # patch contig name in the chr
    patch_gff_contig(
        chr_map,
        'GCF_000149245.1_CNA3_genomic.gff',
        'GCF_000149245.1_CNA3_genomic.patched.gff')
