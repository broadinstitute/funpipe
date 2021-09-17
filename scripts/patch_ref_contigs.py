#!/usr/bin/env python3
import re
import gzip
import sys
import argparse
from funpipe.utils import cd
from funpipe.fasta import samtools_index_fa

'''
    Patch contig names in GFF and fasta from NCBI
    author: Xiao Li (xiaoli@broadinstitute.org)
'''


def chr_map_from_desc(desc_line, ctg_sufx):
    ''' extract contig and chromosome name from description line '''
    if 'mitochondrion' in desc_line:
        is_mt = True
        matched = re.match(r'>(\S+) .+', desc_line)
        contig = matched.group(1)
        chr = 'chrMT'
    elif 'chromosome':
        is_mt = False
        matched = re.match(r'>(\S+) .+ chromosome (\d+),.+', desc_line)
        contig, chr = matched.group(1), matched.group(2)
        if not chr:
            raise ValueError("contig name not matched!")
        chr = 'chr'+chr
    else:
        raise ValueError("Unrecognized description line:" + desc_line)
    chr += ctg_sufx  # add suffix to contig name
    return (contig, chr)


def patch_fasta(in_fa, out_fa, ctg_sufx, rm_mt=True):
    '''
    :param in_fa: input fasta file, can be gziped with suffix gz
    :param out_fa: output fasta file
    :param ctg_sufx: contig suffix
    :param rm_mt: boolean, whether to remove mitochondrion sequence or not
    '''
    chr_map = {}
    is_mt = False
    if in_fa[-2:] == 'gz':
        fa = gzip.open(in_fa, 'rt')
    else:
        fa = open(in_fa, 'r')
    with open(out_fa, 'w') as pfa:
        for line in fa:
            if line[0] == '>':
                (contig, chr) = chr_map_from_desc(line, ctg_sufx)
                chr_map[contig] = chr
                line = '>'+chr+'\n'
                if chr == 'chrMT'+ctg_sufx:
                    is_mt = True
                else:
                    is_mt = False
            if rm_mt:
                if not is_mt:
                    pfa.write(line)
            else:
                pfa.write(line)
    fa.close()
    return chr_map


def patch_gff_contig(chr_map, in_gff, out_gff, ctg_sufx, rm_mt=True):
    '''
    :param chr_map: a dictionary that maps
    :param in_gff: input gff file
    :param out_gff: output gff file
    :param rm_mt: boolean, whether remove mitochondrion genes or not
    '''
    is_mt = False
    with open(in_gff, 'r') as gff, open(out_gff, 'w') as pgff:
        for line in gff:
            if line[0] == '#':   # fix annotations
                if line[:5] == '##seq':
                    seq_region = line.split(' ')
                    seq_region[1] = chr_map[seq_region[1]]
                    if seq_region[1] == 'chrMT'+ctg_sufx:
                        is_mt = True
                    else:
                        is_mt = False
                    line = ' '.join(seq_region)
            else:  # fix genomic features
                gff_line = line.split('\t')
                gff_line[0] = chr_map[gff_line[0]]
                if gff_line[0] == 'chrMT'+ctg_sufx:
                    is_mt = True
                else:
                    is_mt = False
                line = '\t'.join(gff_line)
            if rm_mt:  # output mt
                if not is_mt:
                    pgff.write(line)
            else:
                pgff.write(line)


def write_chr_map(chr_map, out_tsv):
    with open(out_tsv, 'w') as out:
        for (k, v) in chr_map.items():
            out.write(k+'\t'+v+'\n')


def run_patch(dir, prefix, ctg_sufx, idx):
    """
    :param dir: working Directory
    :param prefix: input file prefix
    :param ctg_sufx: contig suffix
    :param idx: whether to index patched fasta file or not.
    """
    with cd(dir):
        chr_map = patch_fasta(
            prefix+'.fna.gz', prefix+'.patched.noMT.fasta', ctg_sufx,
            rm_mt=True)
        if idx:
            samtools_index_fa(prefix+'.patched.noMT.fasta')
        write_chr_map(chr_map, prefix+'_chr_map.tsv')
        # patch contig name in the chr
        patch_gff_contig(
            chr_map, prefix+'.gff', prefix+'.patched_noMT.gff', ctg_sufx)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Patch reference genomes fasta and gff files')

    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-d', '--dir', default='.', help='Working Directory',
        required=True)
    required.add_argument(
        '-p', '--prefix', help='prefix for input and output files',
        required=True)

    # optional arguments
    parser.add_argument(
        '-c', '--ctg_sufx', default='',
        help="Suffix of contig names to distinguish different subgenomes")

    args = parser.parse_args()

    run_patch(args.dir, args.prefix, args.ctg_sufx, idx=True)
