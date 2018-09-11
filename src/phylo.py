from .utils import run


def vcf_snp_to_fasta(invcf, prefix, max_amb=100000):
    ''' snp only vcf to fasta file
    :param invcf: input vcf file
    :param prefix: output file prefix
    :param max_amb: maximum number of samples with ambiguous calls for a site
                    to be included, recommended number of samples 10%
    '''
    cmd = ' '.join(['vcfSnpsToFasta.py --max_amb_samples', max_amb, invcf, '>',
                   prefix+'.fasta'])
    run(cmd)
    return prefix+'.fasta'


def fa2phylip(fa, out_prefix, jar):
    ''' transfer fasta file to phylip with java tool readSeq
    :param fa: fasta file
    :param jar: path to readseq.jar
    :param out_prefix:
    '''
    cmd = ' '.join(['java -cp', jar, 'run -f 12', fa])
    run(cmd)
    return


def ramxl(phylip, output, threads):
    ''' Run RAaML
    :param phylip: input phylip format file
    :param output: output file name
    :param threads: number of threads used for
    '''
    cmd = ' '.join([
        'raxmlHPC-PTHREADS-SSE3 -p 78960 -f a -x 12345 -N 1000 -m GTRCAT',
        '-T', str(threads), '-n', output, '-s', phylip])
    run(cmd)
    return output


def fasttree(fa, prefix):
    ''' Run FastTreeDP
    :param fa: fasta file
    :param prefix: output prefix
    '''
    cmd = ' '.join([
        'FastTreeDP -nt', fa, '>', out_prefix+'.nwk'
    ])
    run(cmd)
    return prefix+'.nwk'


def FastTreeDP(in_fa, out_prefix):
    ''' perform fastaTreeDP analysis
    :param in_fa: input fasta file
    :param out_prefix: output file prefix
    :returns nwk file
    '''
    out_nwk = out_prefix+'.nwk'
    cmd = 'FastTreeDP -nt '+in_fa+' > ' + out_nwk
    run(cmd)
    return out_nwk
